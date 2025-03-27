# FluQC: Long-Read Influenza Quality Control

FluQC is a QC platform specifically designed for influenza sequencing data generated by long-read sequencers, such as the Nanopore platform.

FluQC generates four different pages:
- Summary
A summary table showing general statistics such as number of reads and bases, percentage of reads mapped to reference, N50 and the minimum depth of coverage.
- Mapping Statistics
Heatmaps for mapping statistics per sample, per segment.
- In-depth Sample view
Bivariate Read length v.s. quality plots, read lengths per segment and depth of coverage histograms per sample.
- DIP's
A (rudimentary) percentage of putative defective interfering particles per segment, per sample. The DIP percentage is defined as the fraction of reads with lengths greater than half the length of segment they're mapped to.

## Dependencies
FluQC is written in python with some external libraries and dependencies:
- [Minimap2](https://github.com/lh3/minimap2)
- [samtools](https://github.com/samtools/samtools)
- [seqtk](https://github.com/lh3/seqtk)
- [biopython](https://biopython.org)
- [dash & plotly](https://dash.plotly.com)
- [pandas](https://pandas.pydata.org)


# Installation

To install FluQC, first clone the repository and create the conda environment:
```bash
git clone https://github.com/Avans-ATLS/fluQC.git
conda env create -n fluqc -f fluQC/env.yaml
```

Then, activate the environment and install miniflu from the wheel in `dist/`.
```bash
conda activate fluqc
python -m pip install fluQC/dist/fluqc-0.3.0-py3-none-any.whl

# Editable installation (for devs)
python -m pip install -e fluQC/
```


# Running the dashboard

FluQC consists of two modules: `preprocessing` and `dashboard`.

To generate the QC results, first run the `preprocessing` module on your data:
```bash
FluQC preprocessing /path/to/fastq/ /path/to/database.fasta /path/to/output/
```

After successful completion, the dashboard can be run using the `dashboard_data.pkl` file:
```bash
FluQC dashboard /path/to/output/dashboard_data.pkl
```

This command starts an interactive dashboard which you can open in a local browser at http://127.0.0.1:8050/.


# Usage
```bash
usage: FluQC [-h] {preprocess,dashboard} ...

Launch a QC dashboard for an influenza sequencing run

options:
  -h, --help            show this help message and exit

commands:
  {preprocess,dashboard}
                        Valid Subcommands
    preprocess          Analyze fastq files and prepare data for dashboard
    dashboard           Launch dashboard

Developed by Sander Boden @ Avans-ATLS (s.boden1@avans.nl)
```

## preprocess
```bash
usage: FluQC preprocess [-h] [--threads THREADS] fastq database outdir

positional arguments:
  fastq              path to directory of fastqs to analyze
  database           path to IRMA database
  outdir             Path to directory to place output

options:
  -h, --help         show this help message and exit
  --threads THREADS  Number of threads
```

## dashboard
```bash
usage: FluQC dashboard [-h] datapath

positional arguments:
  datapath    path to directory containing data csv's. (outdir from preprocess subcommand)

options:
  -h, --help  show this help message and exit
```



