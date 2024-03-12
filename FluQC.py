import argparse
import glob
import os
import logging.config
import yaml


import pandas as pd
from dash import Dash, html, dcc, callback, Output, Input
import plotly.express as px
from plotly.graph_objects import Figure, Scatter

from analysis import Analysis

class Dashboard():
    segments = ['HA', 'NA', 'PB1', 'PB2', 'MP', 'NP', 'PA', 'NS']
    
    def plot_percent_dip(self):
        d = dict(
        sample = [],
        HA=[], NA=[], PB1=[], PB2=[], MP=[], NP=[], PA=[], NS=[]
        )
        # combine results into dataframe
        for sample, data in results.items():
            d['sample'].append(sample)
            for k, v in data.percent_dip.items():
                d[k.split('_')[1]].append(v)
             # if length of a key in d is less than len(d[sample]), append 0
            for key in list(d.keys())[1:]:
                if len(d[key]) < len(d['sample']):
                    d[key].append(0)
        df = pd.DataFrame(d)
        df = df.set_index('sample')
        print(df)
        return px.imshow(df, text_auto=True)
    
    
    @callback(Output("heatmap", "figure"), Input("statistic", "value"))
    def create_depth_heatmap_data(value):
        d = dict(
        sample = [],
        HA=[], NA=[], PB1=[], PB2=[], MP=[], NP=[], PA=[], NS=[]
        )
        # iterate over all analyzed samples
        for sample, data in results.items():
            d['sample'].append(sample)
            # for each segment that passed QC checks in analysis:
            for segment in data.segments:
                # add data to df
                d[segment.split('_')[1]].append(data.covstats.loc[segment, value])
            # if length of a key in d is less than len(d[sample]), append 0
            for key in list(d.keys())[1:]:
                if len(d[key]) < len(d['sample']):
                    d[key].append(0)
        df = pd.DataFrame(d)
        df = df.set_index('sample')
        print(df)
        return px.imshow(df, text_auto=True)

    @callback(Output("segment-lengths", "figure"), Input("sample", "value"))
    def plot_segment_lengths(value: str) -> Figure:
        fig =  px.violin(results[value].segment_readlengths, x="segment", y="length", color='segment')

        fig.add_trace(Scatter(
            x=list(Analysis.segment_lengths.keys()), y=list(Analysis.segment_lengths.values()),
            mode='markers', fillcolor='red', name='Expected segment size'                  
                              ))
        return fig

    def __init__(self):
        l = logging.getLogger("Dashboard")
        l.info("Starting Dashboard")
        app = Dash(__name__)
        app.layout = html.Div(
        [
            html.H1(children="FluQC dashboard", style={"textAlign": "center"}),
            html.H3(children='Heatmap of percentage putative DIPs'),
            dcc.Graph(figure=self.plot_percent_dip()),
            html.H2(children='Choose sample to show results:'),
            dcc.Dropdown(['coverage', 'numreads', 'covbases', 'meanbaseq', 'meanmapq', 'meandepth'], 'meandepth', id='statistic'),
            html.H3(children="Heatmap of samtools coverage stats"),
            dcc.Graph(id='heatmap'),
            dcc.Dropdown(list(results.keys()), list(results.keys())[0], id="sample"),
            html.H3(children="Violin plot of segment lengths") ,
            dcc.Graph(id="segment-lengths"),
        ]
    )
        # self.create_depth_heatmap_data()
        # self.plot_percent_dip()
        app.run(debug=True)



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
    paths = [x for x in glob.glob(os.path.join(args.fastq, "*.fastq"))]
    # get a list of samplenames
    names = [os.path.basename(x).replace(".fastq", "") for x in paths]

    # analyze samples
    results = {}
    for p, n in zip(paths, names):
        results[n] = Analysis(p, args.database, args.threads, args.outdir)

    Dashboard()


