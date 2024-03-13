import argparse
import glob
import os
import logging.config
import yaml


import pandas as pd
from pandas import DataFrame
from dash import Dash, html, dcc, callback, Output, Input
import plotly.express as px
from plotly.graph_objects import Figure, Scatter
import plotly.graph_objects as go

from analysis import Analysis


segment_lengths: dict = dict(
        HA=1700, NA=1413, PB1=2274, PB2=2280, MP=982, NP=1497, PA=2151, NS=863
    )



def plot_percent_dip():
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
    df = df.transpose() 
    fig = px.imshow(df, text_auto=True)
    fig.update_xaxes(categoryorder='category ascending')
    return fig

@callback(Output("depth", "figure"), Input("sample", "value"))
def plot_depth(value):
    data = results[value].depth
    
    # aggregate data per 50 bp
    histdf = dict(
        segment = [],
        position = [],
        avg_depth = [],
    )
    for segment in data.segment.unique():
        subdf = data[data['segment'] == segment]
        rolling_avg: pd.Series = subdf['depth'].rolling(window=50).mean().iloc[::50].dropna().reset_index(drop=True)
        for avg in rolling_avg.items():
            histdf["segment"].append(segment.split('_')[1])
            histdf['position'].append((avg[0]+1)*50)
            histdf["avg_depth"].append(avg[1])
    df = pd.DataFrame(histdf)
    
    # make array for multicategory axis
    x = [
        [x for x in df.segment],
        [x for x in df.position],
    ]
    fig = go.Figure()
    fig.add_bar(x=x, y=df.avg_depth, opacity=1, marker_color='MediumPurple')
    fig.update_xaxes(categoryorder='category ascending')
    return fig

@callback(Output("linedepth", "figure"), Input("sample", "value"))    
def plot_depth_line(value):
    data = results[value].depth
    
    # multicategory axis
    x = [
        [x.split('_')[1] for x in data.segment],
        [x for x in data.pos]
    ]
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(x=x, y=data.depth)
    )
    return fig
    
    
@callback(Output("heatmap", "figure"), Input("statistic", "value"))
def plot_samstats_heatmap(value):
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
    df = df.transpose()
    fig = px.imshow(df, text_auto=True)
    fig.update_xaxes(categoryorder='category ascending')
    return fig

@callback(Output("segment-lengths", "figure"), Input("sample", "value"))
def plot_segment_lengths(value: str) -> Figure:
    fig =  px.violin(results[value].segment_readlengths, x="segment", y="length", color='segment')

    fig.add_trace(Scatter(
        x=list(Analysis.segment_lengths.keys()), y=list(Analysis.segment_lengths.values()),
        mode='markers', fillcolor='red', name='Expected segment size'                  
                            ))
    return fig




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

    logger.info(f"Starting analysis of samples in {args.fastq}")
    logger.debug(f'samples: {'\t'.join(names)}')
    # analyze samples
    results = {}
    for p, n in zip(paths, names):
        results[n] = Analysis(p, args.database, args.threads, args.outdir)

    logger = logging.getLogger("Dashboard")
    logger.info("Starting Dashboard")
    app = Dash(__name__)
    app.layout = html.Div(
        [
            html.H1(children="FluQC dashboard", style={"textAlign": "center"}),
            html.H3(children='Heatmap of percentage putative DIPs'),
            dcc.Graph(figure=plot_percent_dip()),
            html.H2(children='Choose sample to show results:'),
            dcc.Dropdown(['coverage', 'numreads', 'covbases', 'meanbaseq', 'meanmapq', 'meandepth'], 'meandepth', id='statistic'),
            html.H3(children="Heatmap of samtools coverage stats"),
            dcc.Graph(id='heatmap'),
            dcc.Dropdown(list(results.keys()), list(results.keys())[0], id="sample"),
            html.H3(children="Violin plot of segment lengths") ,
            dcc.Graph(id="segment-lengths"),
            html.H3(children="Read depth across all segments"),
            dcc.Graph(id="depth"),
            dcc.Graph(id="linedepth")
        ]
    )
    app.run(debug=True)


