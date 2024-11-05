import logging

from dash import Dash, html, dcc, callback, Output, Input
import dash_bootstrap_components as dbc
from plotly.graph_objects import Figure
import pandas as pd
from pandas import DataFrame

from fluqc.figuredata import FigureData
from fluqc.plotter import Plots
from fluqc.dashboard.texts import DashboardText


def read_datatable(d: dict) -> DataFrame:
    """Read

    Args:
        d (dict): dictionary of data (from Figuredata)

    Returns:
        DataFrame: data to show in dashboard table
    """
    # create df, transpose to set samples to cols
    table = pd.DataFrame(d, columns=None)
    # table: pd.DataFrame = table.transpose()
    # # set column labels as sample, drop sample row
    # table.columns = table.loc["Sample"].values
    # table.drop("Sample", inplace=True)
    # table.reset_index()
    # # create a new column from the df index
    # table["Statistic"] = table.index
    # cols = list(table.columns)
    # # move the statistics column from last to first
    # table = table[[cols[-1]] + cols[:-1]]
    # # transpose the table
    # table = table.transpose()

    return table


SIDEBAR = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
}
CONTENT = {
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
}


def launch_dashboard(data: FigureData) -> None:
    """Create html layout and launch QC dashboard

    Args:
        data (FigureData): instance of FigureData class containing data to plot
    """
    app = Dash(
        "FluQC",
        external_stylesheets=[dbc.themes.PULSE],
        suppress_callback_exceptions=True,
    )  # spacelab/pulse are nice

    # callbacks for interaction
    @callback(Output("samstats", "figure"), [Input("statistic", "value")])
    def show_covstats(value) -> Figure:
        return p.cov[value]

    @callback(Output("lengths", "figure"), Input("sample", "value"))
    def show_lenghts(value) -> Figure:
        return p.lengths[value]

    @callback(Output("depth", "figure"), Input("sample", "value"))
    def show_depth(value) -> Figure:
        return p.depth[value]

    @callback(Output("bivariate", "figure"), Input("sample", "value"))
    def show_bivariate(value) -> Figure:
        return p.bivariate[value]

    p = Plots(data)
    t = DashboardText()

    logger = logging.getLogger("Dashboard")
    logger.info("Starting Dashboard")

    sidebar_element = html.Div(
        [
            html.H1("FluQC", className="text-primary"),
            html.Hr(),
            html.P("Avans-ATLS", className="text-primary-emphasis"),
            dbc.Nav(
                [
                    dbc.NavLink("Home", href="/", active="exact"),
                    dbc.NavLink("Summary", href="/page-summary", active="exact"),
                    dbc.NavLink(
                        "Mapping Statistics", href="/page-mapping", active="exact"
                    ),
                    dbc.NavLink(
                        "In-depth Sample view", href="/page-sample", active="exact"
                    ),
                    dbc.NavLink("DIP's", href="/page-dip", active="exact"),
                ],
                vertical=True,
                pills=True,
            ),
        ],
        style=SIDEBAR,
    )

    HOMEPAGE = html.Div(
        [
            html.Br(),
            html.H1("FluQC dashboard", className="text-primary"),
            html.Br(),
            html.H2("Introduction:"),
            dcc.Markdown(DashboardText.introduction),
            html.Br(),
            html.H4("Epilogue"),
            dcc.Markdown(DashboardText.epilogue),
        ],
        style=CONTENT,
    )
    SUMMARY_PAGE = html.Div(
        [
            html.Br(),
            dcc.Markdown(t.table_text),
            html.Br(),
            dbc.Table.from_dataframe(
                read_datatable(data.table), striped=True, bordered=True, hover=True
            ),
        ],
        style=CONTENT,
    )
    MAPPING_PAGE = html.Div(
        [
            html.H1("Mapping Statistics"),
            html.Br(),
            dcc.Markdown(t.mapping_explanation),
            html.Br(),
            html.H4(children="Choose a statistic to see the results"),
            dcc.Dropdown(list(p.cov.keys()), list(p.cov.keys())[0], id="statistic"),
            html.Br(),
            dcc.Graph(id="samstats"),
        ],
        style=CONTENT,
    )
    SAMPLE_PAGE = html.Div(
        [
            html.Br(),
            html.H1("In-depth Sample view"),
            html.Br(),
            html.H4("Choose a sample to view figures"),
            dcc.Dropdown(
                list(p.bivariate.keys()), list(p.bivariate.keys())[0], id="sample"
            ),
            html.Br(),
            html.H5("Read length v.s average quality"),
            dcc.Graph(id="bivariate"),
            html.Br(),
            html.H5("Read lengths per sample"),
            dcc.Graph(id="lengths"),
            html.Br(),
            html.H5("Read depth histogram"),
            dcc.Graph(id="depth"),
        ],
        style=CONTENT,
    )
    DIP_PAGE = html.Div(
        [
            html.H1("Percentage Defective interfering particles"),
            html.Br(),
            dcc.Markdown(t.dip_explanation),
            html.Br(),
            html.H5("DIP Heatmap"),
            html.Br(),
            dcc.Graph(figure=p.dip),
        ],
        style=CONTENT,
    )

    content = html.Div(id="page-content")
    app.layout = html.Div(  # new version
        [
            dcc.Location(id="url"),
            sidebar_element,
            content,
        ]
    )

    # callback for showing pages
    @app.callback(Output("page-content", "children"), [Input("url", "pathname")])
    def render_page(pathname):
        if pathname == "/":
            return HOMEPAGE
        elif pathname == "/page-summary":
            return SUMMARY_PAGE
        elif pathname == "/page-mapping":
            return MAPPING_PAGE
        elif pathname == "/page-sample":
            return SAMPLE_PAGE
        elif pathname == "/page-dip":
            return DIP_PAGE
        else:
            return html.Div(
                [
                    html.H1("404: Not found", className="text-danger"),
                    html.Hr(),
                    html.P(f"The pathname {pathname} was not recognised..."),
                ],
                style=CONTENT,
            )

    app.run(debug=True)
