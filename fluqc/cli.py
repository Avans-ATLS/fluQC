import os
import glob
import pickle
import argparse
import logging.config

import yaml

from fluqc.wrappers import Wrappers
from fluqc.figuredata import FigureData, FastqStats
from fluqc.samplepaths import SamplePaths
from fluqc.dashboard.layout import launch_dashboard


def write_figuredata(data: FigureData, outfile: str):
    with open(outfile, "wb") as out:
        pickle.dump(data, out, pickle.HIGHEST_PROTOCOL)


def run_preprocessing(args):
    logger = logging.getLogger("FluQC")
    # init SamplePaths classes for each fastq file
    samples = [
        SamplePaths(x, args.database, args.outdir)
        for x in sorted(glob.glob(os.path.join(args.fastq, "*.fastq")))
    ]
    logger.info(f"Starting analysis for {len(samples)} samples in {args.fastq}")

    # init FigureData class
    data = FigureData(args.threads)

    # analyze each dataset and collect results
    for s in samples:
        logger.info(f"Analyzing {s.samplename}")
        # check if sample has > 10k reads.
        if FastqStats(s.fastq).nreads > args.sample and not args.sample == 0:
            logger.info(f"{s.samplename} has > {args.sample} reads")
            # subsample to 10k reads
            s.subsample_fastq(args.sample)

        # init Wrappers class for sample
        run = Wrappers(s, args.threads)
        # run analyses
        run.minimap2()
        run.samtools_cov()
        run.samtools_depth()
        segments, assignments = run.assign_reads()

        # check if assignment was unsuccessful
        if segments == None and assignments == None:
            continue  # skip current sample
        # append results to FigureData class
        data.append_percent_dips(s.samplename, s.paf, segments)
        data.append_covstats(s.samplename, s.samcov, segments)
        data.append_depth(s.samplename, s.samdepth, segments)
        data.append_segment_readlengths(s.samplename, s.fastq, assignments)
        data.append_len_qual(s.samplename, s.fastq)
        data.append_table_data(s.samplename, s.paf, s.fastq, s.samdepth)
    write_figuredata(data, os.path.join(args.outdir, "dashboard_data.pkl"))


def dashboard(args):
    with open(args.datapath, "rb") as input:
        data = pickle.load(input)
    launch_dashboard(data)


def main():
    cfg = os.path.join(os.path.dirname(__file__), "config/logging_config.yml")
    with open(cfg, "rt") as f:
        config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)

    root_parser = argparse.ArgumentParser(
        "FluQC",
        description="Launch a QC dashboard for an influenza sequencing run",
        epilog="Developed by Sander Boden @ Avans-ATLS (s.boden1@avans.nl)",
    )
    subparsers = root_parser.add_subparsers(title="commands", help="Valid Subcommands")
    prep_parser = subparsers.add_parser(
        "preprocess", help="Analyze fastq files and prepare data for dashboard"
    )
    prep_parser.add_argument(
        "fastq", type=str, help="path to directory of fastqs to analyze"
    )
    prep_parser.add_argument("database", type=str, help="path to IRMA database")
    prep_parser.add_argument(
        "outdir", type=str, help="Path to directory to place output"
    )
    prep_parser.add_argument("--threads", type=int, default=2, help="Number of threads")
    prep_parser.add_argument(
        "--sample",
        type=int,
        default=0,
        help="Subsample to n reads if sample has more than this number of reads, 0 for no subsampling",
    )
    prep_parser.set_defaults(func=run_preprocessing)

    dash_parser = subparsers.add_parser("dashboard", help="Launch dashboard")
    dash_parser.add_argument(
        "datapath",
        type=str,
        help="path to dashbaord data file (output from preprocess)",
    )
    dash_parser.set_defaults(func=dashboard)

    args = root_parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
