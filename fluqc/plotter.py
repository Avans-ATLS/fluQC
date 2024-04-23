import logging
from functools import partial
from dash import Dash, html, dcc, callback, Output, Input
import plotly.express as px
from plotly.graph_objects import Figure, Scatter
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pandas import DataFrame

from fluqc.figuredata import FigureData


class Plots:
    """Generate Figures for dashboard"""
    segment_lengths = dict(
        HA=1700, NA=1413, PB1=2274, PB2=2280, MP=982, NP=1497, PA=2151, NS=863
    )
    colors = {
        "background": "#FFF7F0",
        "h1": "#009498",
        "h2": "#F2852B",
        "red": "#C6002A",
        "text": "#000000",
        "white": "#FFFFFF",
    }

    def __init__(self, data: FigureData):
        self.d = data
        self.dip: Figure = self.dip_heatmap()
        self.cov: dict[str, Figure] = self.covstats()
        self.depth: dict[str, Figure] = self.cov_histogram()
        self.lengths: dict[str, Figure] = self.readlengths()

    def dip_heatmap(self) -> Figure:
        """Render heatmap of DIP percentages

        Returns:
            Figure: DIP heatmap figure
        """
        df = DataFrame(self.d.dips).set_index("Sample")
        df = df.transpose()
        for col in df.columns:
            df[col] = df[col].astype(float)
        fig = px.imshow(df, text_auto=True)
        fig.update_layout(
            title=f"Heatmap of % putative DIPs",
            xaxis_title="Sample",
            yaxis_title="Segment",
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
        )
        return fig

    def covstats(self) -> dict[str, Figure]:
        """Render covstats heatmaps per statistic

        Returns:
            dict[str, Figure]: statistic: heatmap figure
        """
        df = DataFrame(self.d.covstats)
        callback_options = df.columns[2:]

        # init dict to map figures to options
        d = {}
        # create figure per option and add to dict
        for opt in callback_options:
            subdf = df[["Sample", "Segment", opt]]
            subdf[opt] = subdf[opt].astype(float)
            subdf = subdf.pivot(index="Segment", columns="Sample", values=opt)
            fig = px.imshow(subdf, text_auto=True)
            fig.update_layout(
                title=f"Heatmap of: {opt}",
                xaxis_title="Sample",
                yaxis_title="Segment",
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
            )
            d[opt] = fig

        return d

    def cov_histogram(self) -> dict[str, Figure]:
        """Render per-segment coverage histogram per sample

        Returns:
            dict[str, Figure]: samplename: coverage histogram
        """
        df = DataFrame(self.d.depth)

        # init dict to map figures to samples
        d = {}
        # create figure per sample and add to dict
        for s in df["Sample"].unique():
            subdf = df[df["Sample"] == s]

            fig = make_subplots(
                rows=1,
                cols=len(df["Segment"].unique()),
                subplot_titles=sorted(subdf["Segment"].unique()),
                shared_xaxes=True,
                shared_yaxes=True,
                vertical_spacing=0,
                horizontal_spacing=0.01,
            )

            for i, segment in enumerate(sorted(subdf["Segment"].unique()), start=1):
                segdf = subdf[subdf["Segment"] == segment]
                segdf["position"] = segdf["Position"].astype(int)
                segdf["RollingAvg"] = segdf["RollingAvg"].astype(float)

                fig.add_trace(
                    go.Bar(
                        x=segdf["Position"],
                        y=segdf["RollingAvg"],
                        name=segment,
                        opacity=1,
                        marker_color=self.colors["h1"],
                    ),
                    row=1,
                    col=i,
                )
                x_range = [0, segdf["Position"].max() + 10]
                y_range = [0, subdf["RollingAvg"].max() + 10]

                fig.update_xaxes(range=x_range, tickangle=45, row=1, col=i)
                fig.update_yaxes(range=y_range, tickangle=45, row=1, col=i)

            fig.update_layout(
                title=f"Segment coverage: {s}",
                xaxis_title="Position",
                yaxis_title="Read depth",
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                font_color=self.colors["text"],
            )
            d[s] = fig

        return d

    def readlengths(self) -> dict[str, Figure]:
        """Render violin plots of readlength per sample

        Returns:
            dict[str, Figure]: samplename: violin plot
        """
        df = DataFrame(self.d.readlengths)

        # init dict to map figures to samples
        d = {}
        # create figure per sample and add to dict
        for s in df["Sample"].unique():
            subdf = df[df["Sample"] == s]
            fig = px.violin(
                subdf,
                x="Segment",
                y="ReadLength",
                color_discrete_sequence=[self.colors["h1"]],
            )
            fig.add_trace(
                Scatter(
                    x=list(self.segment_lengths.keys()),
                    y=list(self.segment_lengths.values()),
                    mode="markers",
                    fillcolor=self.colors["h1"],
                    name="Expected segment size",
                )
            )
            fig.update_layout(
                title=f"Segment depth: {s}",
                xaxis_title="Position",
                yaxis_title="Read depth",
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                font_color=self.colors["text"],
            )
            d[s] = fig

        return d
