import argparse
import glob
import os
import logging.config
import yaml
import pickle

from fluqc.figuredata import FigureData
from fluqc.samplepaths import SamplePaths
from fluqc.wrappers import Wrappers, kmerfreqs, kmerPCA
from fluqc.dashboard.layout import launch_dashboard


def get_subtype(segments: list) -> str:
    """From a list of unique assigned reference segments, assign the subtype.
    If HA and/or NA subtype cannot be assigned, Hx and Nx will be returned, respectively.

    Args:
        segments (list): assigned segments

    Returns:
        str: HxNx, where x is either a subtype number, or x if not assigned
    """
    try:
        ha = [x for x in segments if x.split("_")[1] == "HA"][0].split("_")[-1]
    except:
        ha = "Hx"
    try:
        na = [x for x in segments if x.split("_")[1] == "NA"][0].split("_")[-1]
    except:
        na = "Nx"
    return f"{ha}{na}"


def write_figuredata(data: FigureData, outfile: str):
    with open(outfile, "wb") as out:
        pickle.dump(data, out, pickle.HIGHEST_PROTOCOL)


def run_preprocessing(args):
    # init SamplePaths classes for each fastq file
    samples = [
        SamplePaths(x, args.database, args.outdir)
        for x in glob.glob(os.path.join(args.fastq, "*.fastq"))
    ]
    logger.info(f"Starting analysis for {len(samples)} samples in {args.fastq}")

    # init FigureData class
    data = FigureData(args.threads)

    # analyze each dataset and collect results
    for s in samples:
        # init Wrappers class for sample
        run = Wrappers(s, args.threads)
        # run analyses
        run.minimap2()
        run.samtools_cov()
        run.samtools_depth()
        segments, assignments = run.assign_reads()

        # get subtype
        subtype = get_subtype(segments)

        # check if assignment was unsuccessful
        if segments == None and assignments == None:
            continue  # skip current sample
        # append results to FigureData class
        data.append_percent_dips(s.samplename, s.paf, segments)
        data.append_covstats(s.samplename, s.samcov, segments)
        data.append_depth(s.samplename, s.samdepth, segments)
        data.append_segment_readlengths(s.samplename, s.fastq, assignments)
        data.append_len_qual(s.samplename, s.fastq)
        data.append_table_data(s.samplename, s.paf, subtype, s.fastq)
        logger.info("Generating Kmerfrequencies")
        kmerfreqs(s.fastq, s.kmerfreqs, k=5)
        kmerdf = kmerPCA(s.kmerfreqs)
        data.prepare_kmer_data(kmerdf, assignments)
    write_figuredata(data, os.path.join(args.outdir, "dashboard_data.pkl"))


def dashboard(args):
    with open(args.datapath, "rb") as input:
        data = pickle.load(input)
    launch_dashboard(data)


if __name__ == "__main__":
    with open("config/logging_config.yml", "rt") as f:
        config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)

    logger = logging.getLogger("FluQC")
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
    prep_parser.set_defaults(func=run_preprocessing)

    dash_parser = subparsers.add_parser("dashboard", help="Launch dashboard")
    dash_parser.add_argument(
        "datapath",
        type=str,
        help="path to directory containing data csv's. (outdir from preprocess subcommand)",
    )
    dash_parser.set_defaults(func=dashboard)

    args = root_parser.parse_args()
    args.func(args)
