import logging
from functools import partial
from dash import Dash, html, dcc, callback, Output, Input, dash_table
import plotly.express as px
from plotly.graph_objects import Figure, Scatter
import plotly.graph_objects as go
import pandas as pd
from pandas import DataFrame

from fluqc.figuredata import FigureData
from fluqc.plotter import Plots
from fluqc.texts import DashboardText


def launch_dashboard(data: FigureData) -> None:
    """Create html layout and launch QC dashboard

    Args:
        data (FigureData): instance of FigureData class containing data to plot
    """

    @callback(Output("samstats", "figure"), Input("statistic", "value"))
    def show_covstats(value) -> Figure:
        return p.cov[value]

    @callback(Output("lengths", "figure"), Input("sample", "value"))
    def show_lenghts(value) -> Figure:
        return p.lengths[value]

    @callback(Output("depth", "figure"), Input("sample", "value"))
    def show_depth(value) -> Figure:
        return p.depth[value]
    
    @callback(Output("bivariate", 'figure'), Input("qc_sample", 'value'))
    def show_bivariate(value) -> Figure:
        return p.bivariate[value]

    colors = {
        "background": "#FFF7F0",
        "h1": "#009498",
        "h2": "#F2852B",
        "red": "#C6002A",
        "text": "#000000",
        "white": "#FFFFFF",
    }

    p = Plots(data)
    t = DashboardText()

    # create table for dashboard
    table = pd.DataFrame(data.table, columns=None)
    # table.set_index('Sample')
    table: pd.DataFrame = table.transpose()
    table.columns = table.loc['Sample'].values
    table.drop('Sample', inplace=True)
    table.reset_index()
    table['Statistic'] = table.index
    cols = list(table.columns)
    table = table[[cols[-1]] + cols[:-1]]


    logger = logging.getLogger("Dashboard")
    logger.info("Starting Dashboard")
    app = Dash(__name__)
    app.layout = html.Div(
        style={"BackgroundColor": colors["background"]},
        children=[
            html.H1(
                children="FluQC dashboard",
                style={
                    "textAlign": "center",
                    "color": colors["h1"],
                    "BackgroundColor": colors["background"],
                },
            ),
            dcc.Markdown(
                children=t.introduction,
                style={
                    "textAlign": "center",
                    "color": colors["text"],
                    "BackgroundColor": colors["background"],
                },
            ),
            dcc.Markdown(
                children=t.table_text,
                style={
                    "textAlign": "left",
                    "color": colors["text"],
                    "BackgroundColor": colors["background"],
                },
            ),
            dash_table.DataTable(table.to_dict('records')),
            dcc.Dropdown(list(p.bivariate.keys()), list(p.bivariate.keys())[0], id="qc_sample"),
            dcc.Graph(id='bivariate'),
            dcc.Dropdown(list(p.cov.keys()), list(p.cov.keys())[0], id="statistic"),
            html.H2(
                children="--- Mapping Statistics ---",
                style={
                    "textAlign": "center",
                    "color": colors["h2"],
                    "BackgroundColor": colors["background"],
                },
            ),
            dcc.Markdown(
                children=t.mapping_explanation,
                style={
                    "textAlign": "left",
                    "color": colors["text"],
                    "BackgroundColor": colors["background"],
                },
            ),
            html.H3(children="Heatmap of samtools coverage stats"),
            dcc.Graph(id="samstats"),
            html.H2(children="Choose sample to show results for:"),
            dcc.Dropdown(
                list(p.lengths.keys()), list(p.lengths.keys())[0], id="sample"
            ),
            html.H3(children="Violin plot of segment lengths"),
            dcc.Graph(id="lengths"),
            html.H3(children="Read depth across all segments"),
            dcc.Graph(id="depth"),
            html.H2(
                children="--- Percentage Putative Defective Interfering Particles ---",
                style={
                    "textAlign": "center",
                    "color": colors["h2"],
                    "BackgroundColor": colors["background"],
                },
            ),
            dcc.Markdown(
                children=t.dip_explanation,
                style={
                    "textAlign": "left",
                    "color": colors["text"],
                    "BackgroundColor": colors["background"],
                },
            ),
            dcc.Graph(figure=p.dip),
        ],
    )
    app.run(debug=True)
