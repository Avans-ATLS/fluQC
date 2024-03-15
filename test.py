import argparse
import glob
import os
import logging.config
import yaml


import pandas as pd
from refactoring import Wrappers, SamplePaths


if __name__ == "__main__":
    with open("logging_config.yml", "rt") as f:
        config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)

    logger = logging.getLogger("FluQC")
    parser = argparse.ArgumentParser(
        "FluQC.py",
        description="Launch a QC dashboard for an influenza sequencing run",
        epilog="Developed by Sander Boden (s.boden1@avans.nl)",
    )
    parser.add_argument("fastq", help="path to directory of fastqs to analyze")
    parser.add_argument("database", help="path to IRMA database")
    parser.add_argument("outdir", help="Path to directory to place output")
    parser.add_argument("--threads", type=int, default=2, help="Number of threads")
    args = parser.parse_args()

    # get list of sample paths
    samples = [SamplePaths(x, args.database, args.outdir) for x in glob.glob(os.path.join(args.fastq, "*.fastq"))]

    logger.info(f"Starting analysis of samples in {args.fastq}")

    for s in samples:
        run = Wrappers(s)
        run.minimap2(2)
        run.samtools_cov()
        run.samtools_depth()
        segments, assignments = run.assign_reads(2)
        
    
    


