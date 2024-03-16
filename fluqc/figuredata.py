import logging
from logging import Logger

import pandas as pd
from Bio import SeqIO

pd.options.mode.copy_on_write = True


class FigureData:
    """Compute and combine results for visualization"""

    l: Logger = logging.getLogger("FigureData")
    DIP_RATIO: float = 2.0  # if < alignment length / segment length: putative dip
    DEPTH_BINSIZE: int = 25  # binsize for depth histogram

    def __init__(self) -> None:
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

    def append_percent_dips(
        self, samplename: str, paf_path: str, segments: list
    ) -> None:

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
        self, samplename: str, covstats_path: str, segments: list
    ) -> None:

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

    def append_depth(self, samplename: str, depth_path: str, segments) -> None:
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
                self.depth["Position"].append((avg[0] + 1) * 50)
                self.depth["RollingAvg"].append(avg[1])

            # If segment not in data, append None
            for k, v in self.depth.items():
                if len(v) < len(self.depth["Sample"]):
                    self.dips[k].append(None)

    def append_segment_readlengths(
        self, samplename: str, fastq_path: str, assignments: dict
    ) -> None:
        def lengths_from_fastq(fastq_path: str) -> dict:
            self.l.debug("Counting readlengths")
            with open(fastq_path) as fq:
                d = {record.id: len(record.seq) for record in SeqIO.parse(fq, "fastq")}
            return d

        # rename keys to segment names in assignments
        assignments = {k: v.split("_")[1] for k, v in assignments.items()}

        lengths = lengths_from_fastq(fastq_path)
        # for each assigned read, append length to data
        for k, v in assignments.items():
            self.readlengths["Sample"].append(samplename)
            self.readlengths["Segment"].append(v)
            self.readlengths["ReadLength"].append(lengths[k])
