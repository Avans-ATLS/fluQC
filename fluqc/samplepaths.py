import os
import logging
import subprocess
from logging import Logger


class SamplePaths:
    """Store paths for a sample"""

    l: Logger = logging.getLogger("SamplePaths")

    def make_directories(self) -> None:
        """If directory does not exist, create it"""
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir)
            self.l.debug(f"Created output directory: {self.outdir}")

    def subsample_fastq(self, nreads: int) -> None:
        """Subsample fastq file to nreads

        Args:
            nreads (int): Number of reads to subsample
        """
        self.l.info(f"Subsampling {self.samplename} to {nreads} reads")
        fastq = os.path.join(self.outdir, f"{self.samplename}_subsampled.fastq")
        cmd = f"seqtk sample {self.fastq} {nreads} > {fastq}"
        subprocess.run(cmd, shell=True)
        self.fastq = fastq

    def __init__(self, fastq: str, db: str, outdir: str):
        """Assign pathdata to class variables, create output directories

        Args:
            fastq (str): Path to fastq file
            db (str): Path to database
            outdir (str): Base output directory
        """
        self.l.info(f"Creating instance for {os.path.basename(fastq)}")
        self.fastq: str = os.path.abspath(fastq)
        self.samplename: str = os.path.basename(fastq).replace(".fastq", "")
        self.outdir = os.path.join(os.path.abspath(outdir), "wd/", self.samplename)
        self.db = os.path.abspath(db)
        self.paf: str = os.path.join(self.outdir, f"{self.samplename}.paf")
        self.outsam: str = os.path.join(self.outdir, f"{self.samplename}.sam")
        self.samcov: str = os.path.join(self.outdir, f"{self.samplename}_cov.tsv")
        self.samdepth: str = os.path.join(self.outdir, f"{self.samplename}_depth.tsv")
        self.unmapped: str = os.path.join(
            self.outdir, f"{self.samplename}_unmapped.fastq"
        )
        self.make_directories()
