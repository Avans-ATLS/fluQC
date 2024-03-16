import logging
from functools import partial
from dash import Dash, html, dcc, callback, Output, Input
import plotly.express as px
from plotly.graph_objects import Figure, Scatter
import plotly.graph_objects as go
from pandas import DataFrame

from fluqc.figuredata import FigureData
from fluqc.plotter import Plots
from fluqc.texts import DashboardText




def launch_dashboard(data: FigureData) -> None:
    
    @callback(Output("samstats", "figure"), Input("statistic", "value"))
    def show_covstats(value) -> Figure:
        return p.cov[value]
    
    @callback(Output("lengths", "figure"), Input("sample", "value"))
    def show_lenghts(value) -> Figure:
        return p.lengths[value]

    @callback(Output("depth", "figure"), Input("sample", "value"))
    def show_depth(value) -> Figure:
        return p.depth[value]

    colors = {
        'background': '#FFF7F0',
        'h1': "#009498",
        'h2': '#F2852B',
        'red': '#C6002A',
        'text': '#000000',
        'white': '#FFFFFF',
    }


    p = Plots(data)
    t = DashboardText()
    


    logger = logging.getLogger("Dashboard")
    logger.info("Starting Dashboard")
    app = Dash(__name__)
    app.layout = html.Div(style={'BackgroundColor': colors["background"]}, children=[
            html.H1(
                children="FluQC dashboard", 
                style={
                    "textAlign": "center",
                    'color': colors["h1"],
                    'BackgroundColor': colors["background"],
                }
            ),
            dcc.Markdown(
                children=t.introduction,
                style={
                    'textAlign': 'center',
                    'color': colors['text'],
                    'BackgroundColor': colors["background"],
                },
            ),
            html.H2(
                children="--- Percentage Putative Differential Interfering Particles ---",
                style={
                    'textAlign': 'center',
                    'color': colors["h2"],
                    'BackgroundColor': colors["background"],
                },
            ),
            dcc.Markdown(
                children=t.dip_explanation,
                style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'BackgroundColor': colors["background"],
                }
            ),
            dcc.Graph(figure=p.dip),
            html.H2(
                children='--- Mapping Statistics ---',
                style={
                    'textAlign': 'left',
                    'color': colors['h2'],
                    'BackgroundColor': colors["background"],
                },
            ),
            dcc.Markdown(
                children=t.mapping_explanation,
                style={
                    'textAlign': 'left',
                    'color': colors['text'],
                    'BackgroundColor': colors["background"],
                },
            ),
            dcc.Dropdown(list(p.cov.keys()), list(p.cov.keys())[0], id='statistic'),
            html.H3(children="Heatmap of samtools coverage stats"),
            dcc.Graph(id='samstats'),
            html.H2(children="Choose sample to show results for:"),
            dcc.Dropdown(list(p.lengths.keys()), list(p.lengths.keys())[0], id="sample"),
            html.H3(children="Violin plot of segment lengths") ,
            dcc.Graph(id="lengths"),
            html.H3(children="Read depth across all segments"),
            dcc.Graph(id="depth"),
        ],
        )
    app.run(debug=True)









