import time
import logging
import subprocess
import multiprocessing
from logging import Logger
from typing import Generator
from functools import partial
from itertools import product
from collections import Counter

import pandas as pd
from Bio import SeqIO
from pandas import DataFrame
from sklearn.decomposition import PCA
from sklearn.preprocessing import RobustScaler

from fluqc.samplepaths import SamplePaths

pd.options.mode.copy_on_write = True



class Wrappers:
    """Collection of CommandLine wrappers"""

    l: Logger = logging.getLogger("Wrappers")
    MIN_COVERAGE: int = 40  # n reads required to be assigned to segment
    MIXED_THRESHOLD: float = (
        0.2  # min fraction of second most abundant segment to be considered mixed
    )

    def __init__(self, samplepaths: SamplePaths, threads: int) -> None:
        self.p = samplepaths
        self.t = threads

    def minimap2(self) -> None:
        """Run Minimap2 for a sample
        Mapped reads are sorted and written to a samfile.
        Sam format is converted to paf and written to a file.

        Args:
            threads (int): Number of threads to use in minimap2
        """
        # read filtering
        # nanofilt = f'NanoFilt -l 400 -q 8 {self.p.fastq}'
        # filt = subprocess.Popen(
        #     nanofilt,
        #     stdout=subprocess.PIPE,
        #     stderr=subprocess.PIPE,
        #     shell=True,
        # )
        # run minimap2
        minimap = f"minimap2 -ax map-ont -t {self.t} {self.p.db} {self.p.fastq}"
        self.l.info(f"Running minimap for read classification")
        map = subprocess.Popen(
            minimap,
            # stdin=filt.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
        )
        # sort
        sort = subprocess.Popen(
            f"samtools sort -O SAM -",
            stdin=map.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
        )
        # split stdout to file and stdout
        tee = subprocess.Popen(
            f"tee {self.p.outsam}",
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
            df["t_len"] = df["t_len"].astype(int)
            df["aln_len"] = df["aln_len"].astype(int)
            df.to_csv(self.p.paf, sep="\t", index=False)

    def samtools_cov(self) -> None:
        """Run samtools coverage, write output to file"""
        self.l.info("Running samtools coverage")
        # run samtools coverage
        cov = subprocess.Popen(
            f"samtools coverage {self.p.outsam}",
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
            df.to_csv(self.p.samcov, sep="\t", index=False)
            self.l.debug(f"Wrote samtools coverage results to {self.p.samcov}")

    def samtools_depth(self) -> None:
        """Run samtools depth, write output to file"""
        self.l.info("Running samtools depth")
        # run samtools depth
        depth = subprocess.Popen(
            f"samtools depth {self.p.outsam}",
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
            df = pd.DataFrame(df)
            df["depth"] = df["depth"].astype(int)
            df["pos"] = df["pos"].astype(int)
            df.to_csv(self.p.samdepth, sep="\t", index=False)
            self.l.debug(f"Wrote samtools depth results to {self.p.samdepth}")

    def assign_read(self, df: DataFrame, read: str) -> tuple[str, str]:
        """Assign read to a segment (used in func assign_reads)

        Args:
            df (DataFrame): PAF output
            read (str): read to classify

        Returns:
            tuple[str,str]: read_id, assigned segment
        """
        rows = df[df["q_name"] == read]
        idx_min = rows["h_mean"].idxmin()
        return read, df.loc[idx_min, "t_name"]

    def assign_reads(self) -> tuple[list[str], dict[str, str]]:
        """From minimap2 PAF output, find top hit for each read by harmonic mean.
        For hits to multiple HA/NA subtypes, determine if mixed.
        if not mixed, get the subtype with most hits
        A unique list of top hits is returned for constructing a dynamic reference genome.

        Args:
            df (DataFrame): DataFrame of PAF output

        Returns:
            tuple[list, dict]: Unique assigned segments, per read assignments
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

        def assign_parallel(df: DataFrame) -> dict[str, str]:
            """Parellelize assignment of reads to segments

            Args:
                df (DataFrame): PAF data

            Returns:
                dict[str, str]: read id: assigned segment
            """
            n_processes = min(self.t, len(df["q_name"].unique()))
            try:
                pool = multiprocessing.Pool(processes=n_processes)
            except ValueError as e:
                self.l.warning(
                    f"No reads to assign, skipping sample {self.p.samplename}"
                )
                return None
            # create function with fixed df arg
            assign_partial = partial(self.assign_read, df)

            # map function to unique reads in parallel
            assignments = {
                k: v for k, v in pool.map(assign_partial, df["q_name"].unique())
            }
            pool.close()
            pool.join()

            return assignments

        df = pd.read_csv(self.p.paf, sep="\t")
        self.l.debug("Computing harmonic means")
        df["h_mean"] = harmonic_mean(df["q_end"] - df["q_start"], df["n_match"])
        self.l.debug("Counting reads per segment")
        start_time = time.time()
        read_assignments = assign_parallel(df)
        if read_assignments == None:
            return None, None  # return none if no reads were assigned
        segment_counts = pd.Series(read_assignments).value_counts().to_dict()
        end_time = time.time()
        print(f"{self.p.samplename} took: {end_time-start_time} seconds")
        self.l.info(
            (
                f"Counting reads for {self.p.samplename} took: {end_time-start_time} seconds"
            )
        )
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
            self.l.info(f"{self.p.samplename} is unmixed")
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
