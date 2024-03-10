import subprocess
import os
import logging.config
from logging import Logger
import multiprocessing
from functools import partial

from Bio import SeqIO
import pandas as pd
from pandas import DataFrame


class Analysis:
    MIN_COVERAGE: int = 40  # n reads required to be assigned to segment
    MIXED_THRESHOLD: float = (
        0.2  # min fraction of second most abundant segment to be considered mixed
    )
    HA_MIXED: bool = False
    NA_MIXED: bool = False
    l: Logger = logging.getLogger("Analysis")
    segment_lengths: dict = dict(
        HA=1700, NA=1413, PB1=2274, PB2=2280, MP=982, NP=1497, PA=2151, NS=863
    )

    def __init__(self, fastq_path: str, db_path: str, threads: int, outdir: str):
        self.t: int = threads
        self.fastq: str = os.path.abspath(fastq_path)
        self.samplename: str = os.path.basename(fastq_path).replace(".fastq", "")
        self.outdir: str = os.path.abspath(outdir)
        self.db: str = os.path.abspath(db_path)
        self.outsam: str = os.path.abspath(os.path.join(outdir, f"{self.samplename}.sam"))

        if not os.path.isdir(self.outdir):
            self.l.debug(f"Outdir does not exist, creating...")
            os.makedirs(self.outdir)

        self.l.info(f"Running analysis for {self.samplename}")
        paf = self.minimap2()
        self.segments, self.assignments = self.assign_reads(paf)
        self.type: str = self.typing(self.segments)
        self.segment_readlengths: dict = self.length_per_segment()
        self.covstats: DataFrame = self.samtools_cov()
        self.depth: DataFrame = self.samtools_depth()

        self.l.info(f"Finished analysis of {self.samplename}, subtype: {self.type}")
        self.l.info(f"HA_mixed: {self.HA_MIXED}")
        self.l.info(f"NA_mixed: {self.NA_MIXED}")

    def minimap2(self) -> DataFrame:
        """Map fastq reads against IRMA database
        Output a DataFrame of the PAF results from minimap2

        Returns:
            DataFrame: dataframe of PAF output
        """

        # run minimap2
        minimap = f"minimap2 -ax map-ont -t {self.t} {self.db} {self.fastq}"
        self.l.info(f"Running minimap for read classification")
        map = subprocess.Popen(
            minimap,
            stdout=subprocess.PIPE,
            shell=True,
        )

        # sort
        sort = subprocess.Popen(
            f"samtools sort -O SAM -",
            stdin=map.stdout,
            stdout=subprocess.PIPE,
            shell=True,
        )
        # split stdout to file and stdout
        tee = subprocess.Popen(
            f"tee {self.outsam}",
            stdin=sort.stdout,
            stdout=subprocess.PIPE,
            shell=True,
        )
        # convert to paf
        paf = subprocess.Popen(
            f"paftools.js sam2paf -",
            stdin=tee.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
        )
        # get output, check for errors
        stdout, stderr = paf.communicate()
        if paf.returncode != 0:
            self.l.error(f"Error {stderr.decode('utf-8')}")
        else:
            stdout = stdout.decode("utf-8").split("\n")
            stdout = [x for x in stdout if not x == ""]
            df = dict(
                q_name=[],
                q_len=[],
                q_start=[],
                q_end=[],
                strand=[],
                t_name=[],
                t_len=[],
                t_start=[],
                t_end=[],
                n_match=[],
                aln_len=[],
                mapq=[],
            )
            for line in stdout:
                cols = line.split()
                for i, k in enumerate(df):
                    df[k].append(cols[i])
            df = pd.DataFrame(df)
            df["q_end"] = df["q_end"].astype(int)
            df["q_start"] = df["q_start"].astype(int)
            df["n_match"] = df["n_match"].astype(int)
            return df

    def assign_read(self, df: DataFrame, read: str) -> tuple:
        """Assign read to a segment (used in func assign_reads)

        Args:
            df (DataFrame): PAF output
            read (str): read to classify

        Returns:
            tuple: segment name of best hit
        """
        subdf = df[df["q_name"] == read]
        subdf.sort_values("h_mean")
        return read, subdf.iloc[0]["t_name"]

    def assign_reads(self, df: DataFrame) -> list:
        """From minimap2 PAF output, find top hit for each read by harmonic mean.
        For hits to multiple HA/NA subtypes, determine if mixed.
        if not mixed, get the subtype with most hits
        A unique list of top hits is returned for constructing a dynamic reference genome.

        Args:
            df (DataFrame): DataFrame of PAF output

        Returns:
            list: Unique segments which reads mapped to
        """

        def harmonic_mean(map_len: int, n_match: int) -> float:
            """Calculate harmonic mean from a single row of PAF data

            Args:
                row (Series): Row of a DataFrame

            Returns:
                float: Harmonic mean
            """
            h_mean = 2 * ((map_len * n_match) / (map_len + n_match))
            return h_mean

        def assign_parallel(df: DataFrame) -> dict:
            n_processes = min(self.t, len(df["q_name"].unique()))
            pool = multiprocessing.Pool(processes=n_processes)

            # create function with fixed df arg
            assign_partial = partial(self.assign_read, df)

            # map function to unique reads in parallel
            assignments = {
                k: v for k, v in pool.map(assign_partial, df["q_name"].unique())
            }
            pool.close()
            pool.join()

            return assignments

        self.l.debug("Computing harmonic means")
        df["h_mean"] = harmonic_mean(df["q_end"] - df["q_start"], df["n_match"])
        self.l.debug("Counting reads per segment")
        read_assignments = assign_parallel(df)
        segment_counts = pd.Series(read_assignments).value_counts().to_dict()

        for k, v in segment_counts.items():
            self.l.debug(f"{k}:\t{v}")

        # filter segments on minimum reads
        filtered_assignments = {
            k: v
            for k, v in segment_counts.items()
            if k not in [k for k, v in segment_counts.items() if v < self.MIN_COVERAGE]
        }

        self.l.debug(
            f"Discarding segments with readcounts less than {self.MIN_COVERAGE}"
        )
        # get lists of HA, NA and all segments
        ha = [k for k in filtered_assignments if k.split("_")[1] == "HA"]
        na = [k for k in filtered_assignments if k.split("_")[1] == "NA"]
        all_segments = list(set(filtered_assignments.keys()))

        # if <=1 of HA and NA exist, return segments as is (unmixed).
        if len(ha) <= 1 and len(na) <= 1:
            self.l.info(f"{self.samplename} is unmixed")
            return all_segments, read_assignments
        # if multiple HA and/or NA exist, determine if mixed.
        if len(ha) > 1:
            self.l.info("Multiple HA found: Determining if HA is mixed")
            readcounts = {k: v for k, v in segment_counts.items() if k in ha}
            segments = [k for k in readcounts.keys()]
            fraction = readcounts[segments[1]] / readcounts[segments[0]]
            if fraction > self.MIXED_THRESHOLD:
                self.l.info("HA is mixed, removing segments for further analysis")
                self.HA_MIXED = True
                # remove HA from all_segments
                all_segments = [s for s in all_segments if s not in ha]
            else:
                self.l.info("HA is unmixed, selecting segment with most reads")
                # keep only subtype with most reads (first in list)
                all_segments = [s for s in all_segments if s not in ha]
                all_segments.append(segments[0])

        if len(na) > 1:
            self.l.info("Multiple NA found: Determining if NA is mixed")
            readcounts = {k: v for k, v in segment_counts.items() if k in na}
            segments = [k for k in readcounts.keys()]
            fraction = readcounts[segments[1]] / readcounts[segments[0]]
            if fraction > self.MIXED_THRESHOLD:
                self.l.info("NA is mixed, removing segments for further analysis")
                self.NA_MIXED = True
                # remove NA from all_segments
                all_segments = [s for s in all_segments if s not in na]
            else:
                self.l.info("NA is unmixed, selecting segment with most reads")
                # keep only subtype with most reads (first in list)
                all_segments = [s for s in all_segments if s not in na]
                all_segments.append(segments[0])

        return all_segments, read_assignments

    def typing(self, assignments: list) -> str:
        """Extract influenza type and subtype from read assignments

        Args:
            assignments (list): Unique segment assignments

        Returns:
            str: type and subtype "(A_HxNx / B)"
        """
        self.l.info("Determining subtype")
        typed = list(set([x.split("_")[0] for x in assignments]))[0]
        if not typed == "B":
            subtype = [
                x.split("_")[-1] for x in assignments if x.split("_")[1] in ["HA", "NA"]
            ]
            return f"A_{''.join(sorted(subtype))}"
        else:
            return typed

    def length_per_segment(self) -> DataFrame:
        def count_read_lengths() -> dict:
            """Return a dictionary of read lengths

            Returns:
                dict: {read_id (str): read_length (int)}
            """
            self.l.info("Counting read lengths")
            lengths = {}
            with open(self.fastq) as fq:
                for record in SeqIO.parse(fq, "fastq"):
                    lengths[record.id] = len(record.seq)
            return lengths

        def create_length_df(assignments: dict, lengths: dict) -> DataFrame:
            """create a dataframe of read lengths and segment assignments

            Args:
                assignments (dict): read to segment assignments
                lengths (dict): read lenghts

            Returns:
                DataFrame: {read: [], segment: [], length: []}
            """
            self.l.info("Creating dataframe of readlengths per segment")
            assignments = {k: v.split("_")[1] for k, v in assignments.items()}
            df = dict(read=[], segment=[], length=[])
            for id in assignments:
                df["read"].append(id)
                df["segment"].append(assignments[id])
                df["length"].append(lengths[id])
            return pd.DataFrame(df)

        self.l.info("Generating read lengths per segment plot")
        lengths = count_read_lengths()
        return create_length_df(self.assignments, lengths)

    def samtools_cov(self) -> DataFrame:
        self.l.info("Running samtools coverage")
        # run samtools coverage
        cov = subprocess.Popen(
            f"samtools coverage {self.outsam}",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        # check if process ran succesfully
        stdout, stderr = cov.communicate()
        if cov.returncode != 0:
            self.l.error(f"Error {stderr.decode('utf-8')}")
        else:
            # parse stdout  into dataframe
            stdout = stdout.decode("utf-8").split("\n")
            stdout = [x.split("\t") for x in stdout[1:] if not x == ""]
            df = dict(
                segment=[],
                startpos=[],
                endpos=[],
                numreads=[],
                covbases=[],
                coverage=[],
                meandepth=[],
                meanbaseq=[],
                meanmapq=[],
            )
            for line in stdout:
                for i, k in enumerate(df):
                    df[k].append(line[i])
            df = pd.DataFrame(df)
            df = df.set_index('segment')
            return df

    def samtools_depth(self) -> DataFrame:
        self.l.info("Running samtools depth")
        # run samtools depth
        depth = subprocess.Popen(
            f"samtools depth {self.outsam}",
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        # check if process ran succesfully
        stdout, stderr = depth.communicate()
        if depth.returncode != 0:
            self.l.error(f"Error {stderr.decode('utf-8')}")
        else:
            # parse stdout  into dataframe
            stdout = stdout.decode("utf-8").split("\n")
            stdout = [x.split("\t") for x in stdout if not x == ""]
            df = dict(
                segment=[],
                pos=[],
                depth=[],
            )
            for line in stdout:
                for i, k in enumerate(df):
                    df[k].append(line[i])
            return pd.DataFrame(df)
