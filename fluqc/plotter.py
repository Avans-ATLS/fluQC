import plotly.express as px
from plotly.graph_objects import Figure
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pandas import DataFrame
import pandas as pd

from fluqc.figuredata import FigureData


class Plots:
    """Generate Figures for dashboard"""

    segment_lengths = dict(
        HA=1700, NA=1413, PB1=2274, PB2=2280, MP=982, NP=1497, PA=2151, NS=863
    )
    segment_order = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
    template = "plotly"
    marker_color = Figure().layout["template"]["layout"]["colorway"][3]

    def __init__(self, data: FigureData):
        self.d = data
        self.dip: Figure = self.dip_heatmap()
        # self.kmer: Figure = self.kmerfreq_PCA()
        self.cov: dict[str, Figure] = self.covstats()
        self.depth: dict[str, Figure] = self.cov_histogram()
        self.lengths: dict[str, Figure] = self.readlengths()
        self.bivariate: dict[str, Figure] = self.len_qual()

    def dip_heatmap(self) -> Figure:
        """Render heatmap of DIP percentages

        Returns:
            Figure: DIP heatmap figure
        """
        df = DataFrame(self.d.dips).set_index("Sample")
        # df["Segment"] = pd.Categorical(df["Segment"], categories=self.segment_order, ordered=True)
        df = df.transpose()
        df = df.reindex(self.segment_order)
        for col in df.columns:
            df[col] = df[col].astype(float)
        fig = px.imshow(df, text_auto=True)
        fig.update_layout(
            title=f"Heatmap of % putative DIPs",
            xaxis_title="Sample",
            yaxis_title="Segment",
            template=self.template,
        )
        return fig

    def len_qual(self) -> dict[str, Figure]:
        """Render length v.s. quality bivariate plots per sample

        Returns:
            dict[str, Figure]: bivariate plot
        """
        df = DataFrame(self.d.len_qual)
        callback_options = df["Sample"].unique().tolist()

        d = {}
        for opt in callback_options:
            subdf = df[df["Sample"] == opt]
            fig = px.scatter(subdf, x="length", y="quality", template=self.template)
            fig.update_layout(
                title=f"Sample: {opt}",
                xaxis_title="Read Length",
                yaxis_title="Read Quality",
            )
            fig.update_traces(marker=dict(color=self.marker_color, size=3, opacity=0.8))
            d[opt] = fig

        return d

    def covstats(self) -> dict[str, Figure]:
        """Render covstats heatmaps per statistic

        Returns:
            dict[str, Figure]: statistic: heatmap figure
        """
        df = DataFrame(self.d.covstats)
        df["Segment"] = pd.Categorical(
            df["Segment"], categories=self.segment_order, ordered=True
        )
        callback_options = df.columns[2:]

        # init dict to map figures to options
        d = {}
        # create figure per option and add to dict
        for opt in callback_options:
            subdf = df[["Sample", "Segment", opt]]
            subdf[opt] = subdf[opt].astype(float)
            subdf = subdf.pivot(index="Segment", columns="Sample", values=opt)
            subdf = subdf.reindex(self.segment_order)
            fig = px.imshow(subdf, text_auto=True)
            fig.update_layout(
                title=f"Heatmap of: {opt}",
                xaxis_title="Sample",
                yaxis_title="Segment",
                template=self.template,
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
            available_segments = subdf["Segment"].unique()
            sorted_segments = [
                seg for seg in self.segment_order if seg in available_segments
            ]
            fig = make_subplots(
                rows=1,
                cols=len(df["Segment"].unique()),
                subplot_titles=sorted_segments,
                shared_xaxes=True,
                shared_yaxes=True,
                vertical_spacing=0,
                horizontal_spacing=0.01,
            )

            for i, segment in enumerate(sorted_segments, start=1):
                segdf = subdf[subdf["Segment"] == segment]
                segdf["position"] = segdf["Position"].astype(int)
                segdf["RollingAvg"] = segdf["RollingAvg"].astype(float)

                fig.add_trace(
                    go.Bar(
                        x=segdf["Position"],
                        y=segdf["RollingAvg"],
                        name=segment,
                        opacity=1,
                        marker_color=self.marker_color,
                    ),
                    row=1,
                    col=i,
                )
                x_range = [0, segdf["Position"].max() + 5]
                y_range = [0, subdf["RollingAvg"].max() + 1]

                fig.update_xaxes(range=x_range, tickangle=45, row=1, col=i)
                fig.update_yaxes(range=y_range, tickangle=45, row=1, col=i)

            fig.update_layout(
                title=f"Segment coverage: {s}",
                xaxis_title="Position",
                yaxis_title="Read depth (log10)",
                template=self.template,
                showlegend=False,
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
            available_segments = df["Segment"].unique()
            sorted_segments = [
                seg for seg in self.segment_order if seg in available_segments
            ]
            subdf["Segment"] = pd.Categorical(
                subdf["Segment"], categories=self.segment_order, ordered=True
            )
            fig = px.violin(
                subdf,
                x="Segment",
                y="ReadLength",
                category_orders={"Segment": sorted_segments},
            )
            fig.update_traces(marker=dict(color=self.marker_color, size=3, opacity=0.8))
            fig.add_trace(
                go.Scatter(
                    x=list(self.segment_lengths.keys()),
                    y=list(self.segment_lengths.values()),
                    mode="markers",
                    fillcolor="rgb(255,0,0)",
                    name="Expected segment size",
                )
            )
            fig.update_layout(
                title=f"Segment depth: {s}",
                xaxis_title="Position",
                yaxis_title="Read length",
                template=self.template,
            )
            fig.update_yaxes(range=[0, max(self.segment_lengths.values()) + 2000])
            d[s] = fig

        return d

    ## Deprecated
    # def kmerfreq_PCA(self):
    #     return px.scatter_3d(
    #         self.d.kmerfreq,
    #         x="PC1",
    #         y="PC2",
    #         z="PC3",
    #         color="mapped_to",
    #         template=self.template,
    #     )
