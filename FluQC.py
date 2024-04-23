import argparse
import glob
import os
import logging.config
import yaml

from fluqc.figuredata import FigureData
from fluqc.samplepaths import SamplePaths
from fluqc.wrappers import Wrappers
from fluqc.dash import launch_dashboard

def get_subtype(segments: list) -> str:
    try:
        ha = [x for x in segments if x.split('_')[1] == 'HA'][0].split('_')[-1]
    except:
        ha = 'Hx'
    try:
        na = [x for x in segments if x.split('_')[1] == 'NA'][0].split('_')[-1]
    except:
        na = 'Nx'
    return f'{ha}{na}'



if __name__ == "__main__":
    with open("config/logging_config.yml", "rt") as f:
        config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)

    logger = logging.getLogger("FluQC")
    parser = argparse.ArgumentParser(
        "FluQC.py",
        description="Launch a QC dashboard for an influenza sequencing run",
        epilog="Developed by Sander Boden (s.boden1@avans.nl)",
    )
    parser.add_argument(
        "fastq", type=str, help="path to directory of fastqs to analyze"
    )
    parser.add_argument("database", type=str, help="path to IRMA database")
    parser.add_argument("outdir", type=str, help="Path to directory to place output")
    parser.add_argument("--threads", type=int, default=2, help="Number of threads")
    args = parser.parse_args()

    # init SamplePaths classes for each fastq file
    samples = [
        SamplePaths(x, args.database, args.outdir)
        for x in glob.glob(os.path.join(args.fastq, "*.fastq"))
    ]
    logger.info(f"Starting analysis for {len(samples)} samples in {args.fastq}")

    # init FigureData class
    data = FigureData(args.threads)

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
            continue # skip current sample
        # append results to FigureData class
        data.append_percent_dips(s.samplename, s.paf, segments)
        data.append_covstats(s.samplename, s.samcov, segments)
        data.append_depth(s.samplename, s.samdepth, segments)
        data.append_segment_readlengths(s.samplename, s.fastq, assignments)
        data.append_table_data(s.samplename, s.paf, subtype, s.fastq)


    launch_dashboard(data)
