import logging
from logging import Logger
import io
import multiprocessing as mp

import pandas as pd
from Bio import SeqIO

pd.options.mode.copy_on_write = True


class FastqStats:
    def __init__(self, fastq_path: str) -> None:
        self.fq = fastq_path
        self.nreads, self.nbases = self.count_reads()
        self.length = self.avg_len()
        self.qual = self.avg_qual()
        self.N50 = self.calculate_read_n50()

    def count_reads(self) -> int:
        nreads: int = 0
        nbases: int = 0
        with open(self.fq) as fq:
            for record in SeqIO.parse(fq, 'fastq'):
                nreads += 1
                nbases += len(record.seq)
        return nreads, nbases
    
    def avg_len(self) -> float:
        total_len: int = 0
        with open(self.fq) as fq:
            for record in SeqIO.parse(fq, 'fastq'):
                total_len += len(record.seq)
        return round(total_len / self.nreads, 1)

    def avg_qual(self) -> float:
        total_qual: int = 0
        with open(self.fq) as fq:
            for record in SeqIO.parse(fq, 'fastq'):
                total_qual += sum(record.letter_annotations['phred_quality'])
        return round(total_qual / self.nbases)

    def calculate_read_n50(self) -> int:
        with open(self.fq) as fq:
            read_lengths = [len(r.seq) for r in SeqIO.parse(fq, 'fastq')]

        # Sort the read lengths in descending order
        sorted_lengths = sorted(read_lengths, reverse=True)
        
        # Calculate the total number of bases
        total_bases = sum(sorted_lengths)
        
        # Calculate the cumulative sum of read lengths
        cumulative_sum = 0
        for length in sorted_lengths:
            cumulative_sum += length
            # Check if cumulative sum exceeds half of the total bases
            if cumulative_sum >= total_bases / 2:
                return length



class FigureData:
    """Compute and combine results for visualization"""

    l: Logger = logging.getLogger("FigureData")
    DIP_RATIO: float = 2.0  # if < alignment length / segment length: putative dip
    DEPTH_BINSIZE: int = 25  # binsize for depth histogram

    def __init__(self, threads: int) -> None:
        self.t = threads
        # create empty dicts:
        self.dips: dict[str, list] = {
            "Sample": [],
            "HA": [],
            "NA": [],
            "PB1": [],
            "PB2": [],
            "MP": [],
            "NP": [],
            "PA": [],
            "NS": [],
        }
        self.covstats: dict[str, list] = {
            "Sample": [],
            "Segment": [],
            "numreads": [],
            "covbases": [],
            "coverage": [],
            "meandepth": [],
            "meanbaseq": [],
            "meanmapq": [],
        }
        self.depth: dict[str, list] = {
            "Sample": [],
            "Segment": [],
            "Position": [],
            "RollingAvg": [],
        }
        self.readlengths: dict[str, list] = {
            "Sample": [],
            "Segment": [],
            "ReadLength": [],
        }
        # create table data
        self.table: dict[str, list] = {
            "Sample": [],
            "Subtype": [],
            "reads": [],
            "bases": [],
            "% mapped": [],
            "avg length": [],
            "read N50": [],
            "avg quality": [],
        }
    def append_table_data(
            self, samplename: str, paf_path: str, subtype: str, fastq_path: str
    ) -> None:
        
        paf = pd.read_csv(paf_path, sep='\t')
        stats = FastqStats(fastq_path)

        mapped_reads = len(paf['q_name'].unique())
        perc_mapped = round(((mapped_reads / stats.nreads) * 100), 1)
        self.table['Sample'].append(samplename)
        self.table['Subtype'].append(subtype)
        self.table['reads'].append(stats.nreads)
        self.table['bases'].append(stats.nbases)
        self.table['% mapped'].append(perc_mapped)
        self.table['avg length'].append(stats.length)
        self.table['read N50'].append(stats.N50)
        self.table['avg quality'].append(stats.qual)
        
        
    def append_percent_dips(
        self, samplename: str, paf_path: str, segments: list[str]
    ) -> None:
        """Append percentage of DIPs per segment to results for a sample
        %DIP is estimated by calculating the ratio of segment length / alignment length.
        Ratio's passing the threshold are counted as a putative DIP.
        A threshold of 2 will count an alignment as DIP when it's alignment lenght is
        less than half the segment lenght

        Args:
            samplename (str): Name of current sample
            paf_path (str): Path to paf file
            segments (list[str]): Assigned segments
        """

        paf = pd.read_csv(paf_path, sep="\t")
        # add current samplename to data
        self.dips["Sample"].append(samplename)
        # compute DIP for each segment and append to results
        for segment in segments:
            # subset paf on current segment
            df = paf[paf["t_name"] == segment]
            # grab number of unique reads in subset
            n_reads = len(df["q_name"].unique())
            # compute alignment ratio (length of target / alignment length)
            df["aln_ratio"] = df.apply(
                lambda row: row["t_len"] / row["aln_len"], axis=1
            )
            # grab count reads where alignment ratio > DIP_RATIO
            n_dip = len(df[df["aln_ratio"] > self.DIP_RATIO]["q_name"].unique())
            # calculate percentage
            p_dip = (n_dip / n_reads) * 100
            # append to results
            self.dips[segment.split("_")[1]].append(round(p_dip, 1))

        # If segment not in data, append None
        for k, v in self.dips.items():
            if len(v) < len(self.dips["Sample"]):
                self.dips[k].append(None)

    def append_covstats(
        self, samplename: str, covstats_path: str, segments: list[str]
    ) -> None:
        """Append results from samtools coverage to results for a sample

        Args:
            samplename (str): Name of current sample
            covstats_path (str): Path to samtools coverage results
            segments (list[str]): Assigned segments
        """

        df = pd.read_csv(covstats_path, sep="\t", index_col=0)

        # for each segment, append data
        for segment in segments:
            self.covstats["Sample"].append(samplename)
            self.covstats["Segment"].append(segment.split("_")[1])

            for stat in df.columns[2:]:
                self.covstats[stat].append(df.loc[segment, stat])

            # If segment not in data, append None
            for k, v in self.covstats.items():
                if len(v) < len(self.covstats["Sample"]):
                    self.dips[k].append(None)

    def append_depth(
        self, samplename: str, depth_path: str, segments: list[str]
    ) -> None:
        """Append data from samtools depth to results for a sample

        Args:
            samplename (str): Name of current sample
            depth_path (str): Path to samtools depth results
            segments (list[str]): Assigned segments
        """
        df = pd.read_csv(depth_path, sep="\t")

        # calculate rolling average per DEPTH_BINSIZE
        for segment in segments:
            subdf = df[df["segment"] == segment]
            rolling_avg = (
                subdf["depth"]
                .rolling(window=self.DEPTH_BINSIZE)
                .mean()
                .iloc[:: self.DEPTH_BINSIZE]
                .dropna()
                .reset_index(drop=True)
            )

            # start position 0 at a coverage of 0
            self.depth["Sample"].append(samplename)
            self.depth["Segment"].append(segment.split("_")[1])
            self.depth["Position"].append(0)
            self.depth["RollingAvg"].append(0)
            # add rolling avg for segment
            for avg in rolling_avg.items():
                self.depth["Sample"].append(samplename)
                self.depth["Segment"].append(segment.split("_")[1])
                self.depth["Position"].append((avg[0] + 1) * self.DEPTH_BINSIZE)
                self.depth["RollingAvg"].append(avg[1])

            # If segment not in data, append None
            for k, v in self.depth.items():
                if len(v) < len(self.depth["Sample"]):
                    self.dips[k].append(None)

    def append_segment_readlengths(
        self, samplename: str, fastq_path: str, assignments: dict[str, str]
    ) -> None:
        """Append readlengths per segment to results for a sample

        Args:
            samplename (str): Name of current sample
            fastq_path (str): Path to fastq file
            assignments (dict[str, str]): Assignments of reads to segments
        """

        def lengths_from_fastq(fastq_path: str) -> dict[str, int]:
            """Get readlenghts from a fastq file

            Args:
                fastq_path (str): Path to fastq file

            Returns:
                dict[str, int]: read id: read length
            """
            self.l.debug("Counting readlengths")
            with open(fastq_path) as fq:
                d = {}
                for record in SeqIO.parse(fq, "fastq"):
                    d[record.id] = len(record.seq)
            return d

        # rename keys to segment names in assignments
        assignments = {k: v.split("_")[1] for k, v in assignments.items()}

        lengths = lengths_from_fastq(fastq_path)
        # for each assigned read, append length to data
        for k, v in assignments.items():
            self.readlengths["Sample"].append(samplename)
            self.readlengths["Segment"].append(v)
            self.readlengths["ReadLength"].append(lengths[k])