#!/usr/bin/env -S uvx marimo run --sandbox --headless --port 8989 --host 0.0.0.0 --no-token
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "altair==5.5.0",
#     "anthropic==0.57.1",
#     "duckdb==1.3.2",
#     "marimo",
#     "pafpy==0.2.0",
#     "polars[pyarrow]==1.31.0",
#     "ruff==0.12.3",
#     "sqlglot==27.0.0",
#     "vegafusion==2.0.2",
#     "vl-convert-python==1.8.0",
# ]
# ///

import marimo

__generated_with = "0.14.10"
app = marimo.App(width="full")


@app.cell
def _():
    import marimo as mo
    import polars as pl
    import altair as alt
    from typing import List, Tuple, Dict, Optional
    from pafpy import PafFile
    import io
    return Optional, PafFile, Tuple, alt, io, mo, pl


@app.cell
def _(mo):
    paf_source = mo.ui.file(label="Select .paf", filetypes=[".paf"])
    return (paf_source,)


@app.cell
def _(mo, paf_source):
    mo.vstack([paf_source, mo.md(f"## {paf_source.value[0].name}") if paf_source.value else mo.md(text="## No .paf loaded.")])
    return


@app.cell
def _(mo):
    interactive = mo.ui.checkbox(label="Zoomable")
    return (interactive,)


@app.cell
def _(alt, interactive, io, load_and_visualize_paf, mo, paf_source):
    df, summary, chart = load_and_visualize_paf(
        io.BytesIO(paf_source.value[0].contents if paf_source.value else b""), width=1200, height=1200
    )

    interval = alt.selection_interval()

    chart = mo.ui.altair_chart(
        chart.interactive() if interactive.value else chart.add_params(interval)
    )

    mo.vstack([interactive, chart, interval]) if paf_source.value else None
    return (chart,)


@app.cell
def _(chart, mo, paf_source):
    table = mo.ui.table(chart.value) if paf_source.value else None
    table
    return (table,)


@app.cell
def _(paf_source, table):
    table.value if paf_source.value else None
    return


@app.cell
def _(Optional, PafFile, Tuple, alt, pl):
    def parse_paf_to_polars(paf_filename: str) -> pl.DataFrame:
        """
        Parse a PAF file using pafpy and convert to a Polars DataFrame

        Args:
            paf_filename: Path to the PAF file

        Returns:
            Polars DataFrame with alignment data
        """
        records = []

        with PafFile(paf_filename) as paf:
            for alignment_id, record in enumerate(paf):
                # Extract basic alignment information
                record_data = {
                    "alignment_id": alignment_id,
                    "query_name": record.qname,
                    "query_length": record.qlen,
                    "query_start": record.qstart,
                    "query_end": record.qend,
                    "target_name": record.tname,
                    "target_length": record.tlen,
                    "target_start": record.tstart,
                    "target_end": record.tend,
                    "strand": str(record.strand),
                    "is_reverse": str(record.strand)
                    == reversed(str(record.strand)),
                    "is_primary": record.is_primary(),
                    "is_secondary": record.is_secondary(),
                    "mapping_quality": record.mapq,
                    "match_length": record.mlen,
                    "block_length": record.blen,
                    "alignment_length": record.qend - record.qstart,
                    "target_alignment_length": record.tend - record.tstart,
                    "query_coverage": record.query_coverage,
                    "target_coverage": record.target_coverage,
                    "blast_identity": record.blast_identity(),
                }

                # Extract common tags if they exist
                for tag_name in [
                    "NM",
                    "AS",
                    "ms",
                    "nn",
                    "tp",
                    "cm",
                    "s1",
                    "s2",
                    "de",
                    "rl",
                ]:
                    try:
                        tag = record.get_tag(tag_name)
                        if tag:
                            record_data[f"tag_{tag_name}"] = tag.value
                    except:
                        record_data[f"tag_{tag_name}"] = None

                records.append(record_data)

        return pl.DataFrame(records)


    def calculate_global_coordinates(df: pl.DataFrame) -> pl.DataFrame:
        """
        Calculate global coordinates for sequences ordered by length
        Handles reverse strand alignments by swapping coordinates

        Args:
            df: DataFrame with alignment data

        Returns:
            DataFrame with added global coordinate columns
        """
        # Get unique sequences and their lengths
        target_info = (
            df.select(["target_name", "target_length"])
            .unique()
            .sort("target_length", descending=True)
        )
        query_info = (
            df.select(["query_name", "query_length"])
            .unique()
            .sort("query_length", descending=True)
        )

        # Calculate cumulative offsets
        target_offsets = target_info.with_columns(
            [
                pl.col("target_length")
                .cum_sum()
                .shift(1, fill_value=0)
                .alias("target_offset")
            ]
        )

        query_offsets = query_info.with_columns(
            [
                pl.col("query_length")
                .cum_sum()
                .shift(1, fill_value=0)
                .alias("query_offset")
            ]
        )

        # Join back to main dataframe
        df_with_coords = df.join(
            target_offsets, on="target_name", how="left"
        ).join(query_offsets, on="query_name", how="left")

        # Calculate global coordinates
        # For reverse strand alignments, swap the query coordinates to draw from right to left
        df_with_coords = df_with_coords.with_columns(
            [
                (pl.col("target_offset") + pl.col("target_start")).alias(
                    "global_target_start"
                ),
                (pl.col("target_offset") + pl.col("target_end")).alias(
                    "global_target_end"
                ),
                # For forward strand: normal coordinates
                # For reverse strand: swap start and end to draw from right to left
                pl.when(pl.col("strand") == "-")
                .then(pl.col("query_offset") + pl.col("query_end"))
                .otherwise(pl.col("query_offset") + pl.col("query_start"))
                .alias("global_query_start"),
                pl.when(pl.col("strand") == "-")
                .then(pl.col("query_offset") + pl.col("query_start"))
                .otherwise(pl.col("query_offset") + pl.col("query_end"))
                .alias("global_query_end"),
                # Add a visual indicator for strand orientation
                pl.when(pl.col("strand") == "-")
                .then(pl.lit("reverse"))
                .otherwise(pl.lit("forward"))
                .alias("strand_orientation"),
            ]
        )

        return df_with_coords


    def get_sequence_boundaries(df: pl.DataFrame) -> pl.DataFrame:
        """
        Get sequence boundary positions for drawing grid lines

        Args:
            df: DataFrame with global coordinates

        Returns:
            DataFrame with boundary information
        """
        # Target boundaries (vertical lines)
        target_bounds = (
            df.select(["target_name", "target_offset"])
            .unique()
            .filter(pl.col("target_offset") > 0)
        )
        target_bounds = target_bounds.with_columns(
            [
                pl.lit("target").alias("boundary_type"),
                pl.lit("vertical").alias("line_type"),
                pl.col("target_offset").alias("position"),
            ]
        ).select(["boundary_type", "line_type", "target_name", "position"])

        # Query boundaries (horizontal lines)
        query_bounds = (
            df.select(["query_name", "query_offset"])
            .unique()
            .filter(pl.col("query_offset") > 0)
        )
        query_bounds = (
            query_bounds.with_columns(
                [
                    pl.lit("query").alias("boundary_type"),
                    pl.lit("horizontal").alias("line_type"),
                    pl.col("query_offset").alias("position"),
                ]
            )
            .select(["boundary_type", "line_type", "query_name", "position"])
            .rename({"query_name": "target_name"})
        )

        # Combine boundaries
        boundaries = pl.concat([target_bounds, query_bounds])

        return boundaries


    def create_dotplot(
        df: pl.DataFrame,
        boundaries_df: Optional[pl.DataFrame] = None,
        width: int = 800,
        height: int = 600,
        title: str = "PAF Alignment Dotplot",
    ) -> alt.Chart:
        """
        Create an interactive dotplot using Vega-Altair
        """
        if df.is_empty():
            return alt.Chart().mark_text(text="No data to display", fontSize=20)

        # Prepare custom tick values for grid alignment
        target_ticks = [0]
        query_ticks = [0]

        if boundaries_df is not None and not boundaries_df.is_empty():
            target_bounds = boundaries_df.filter(pl.col("line_type") == "vertical")
            query_bounds = boundaries_df.filter(
                pl.col("line_type") == "horizontal"
            )

            if not target_bounds.is_empty():
                target_ticks.extend(target_bounds["position"].to_list())
            if not query_bounds.is_empty():
                query_ticks.extend(query_bounds["position"].to_list())

        # Use mark_rule instead of mark_line to avoid double legend
        chart = (
            alt.Chart(df)
            .mark_rule(strokeWidth=3, opacity=0.7)
            .encode(
                x=alt.X(
                    "global_target_start:Q",
                    title="Target Position (Global)",
                    scale=alt.Scale(nice=False),
                    axis=alt.Axis(values=target_ticks)
                    if target_ticks
                    else alt.Axis(),
                ),
                y=alt.Y(
                    "global_query_start:Q",
                    title="Query Position (Global)",
                    scale=alt.Scale(nice=False),
                    axis=alt.Axis(values=query_ticks)
                    if query_ticks
                    else alt.Axis(),
                ),
                x2=alt.X2("global_target_end:Q"),
                y2=alt.Y2("global_query_end:Q"),
                color=alt.Color(
                    "query_name:N",
                    title="Query Sequence",
                    scale=alt.Scale(scheme="category20"),
                ),
                opacity=alt.condition(
                    alt.datum.is_primary, alt.value(0.8), alt.value(0.4)
                ),
                tooltip=[
                    "query_name:N",
                    "target_name:N",
                    "strand:N",
                    "strand_orientation:N",
                    "global_query_start:Q",
                    "global_query_end:Q",
                    "global_target_start:Q",
                    "global_target_end:Q",
                    "alignment_length:Q",
                    "blast_identity:Q",
                    "query_coverage:Q",
                    "mapping_quality:Q",
                    "is_primary:N",
                ],
            )
        )

        # Configure the chart
        chart = chart.resolve_scale(color="independent").properties(
            width=width,
            height=height,
            title=alt.TitleParams(text=title, fontSize=16, fontWeight="bold"),
        )

        return chart


    def create_alignment_summary(df: pl.DataFrame) -> pl.DataFrame:
        """
        Create a summary of alignment statistics

        Args:
            df: DataFrame with alignment data

        Returns:
            Summary DataFrame
        """
        summary = (
            df.group_by(["query_name", "target_name"])
            .agg(
                [
                    pl.len().alias("num_alignments"),
                    pl.col("blast_identity").mean().alias("avg_identity"),
                    pl.col("blast_identity").max().alias("max_identity"),
                    pl.col("query_coverage").sum().alias("total_query_coverage"),
                    pl.col("alignment_length")
                    .sum()
                    .alias("total_alignment_length"),
                    pl.col("is_primary").sum().alias("num_primary"),
                ]
            )
            .sort("total_alignment_length", descending=True)
        )

        return summary


    def load_and_visualize_paf(
        paf_filename: str, width: int = 800, height: int = 600
    ) -> Tuple[pl.DataFrame, pl.DataFrame, alt.Chart]:
        """
        Main function to load PAF file and create visualization using pafpy

        Args:
            paf_filename: Path to PAF file
            width: Chart width in pixels
            height: Chart height in pixels

        Returns:
            Tuple of (alignment DataFrame, summary DataFrame, Altair chart object)
        """

        # Parse PAF file using pafpy
        df = parse_paf_to_polars(paf_filename)

        if df.is_empty():
            print("No alignment records found!")
            return (
                df,
                pl.DataFrame(),
                alt.Chart().mark_text(text="No data to display"),
            )

        # Calculate global coordinates
        df_with_coords = calculate_global_coordinates(df)

        # Get sequence boundaries
        boundaries_df = get_sequence_boundaries(df_with_coords)

        # Create summary
        summary_df = create_alignment_summary(df_with_coords)

        # Create visualization
        chart = create_dotplot(
            df_with_coords,
            boundaries_df,
            width=width,
            height=height,
            title=f"PAF Alignment Dotplot",
        )

        return df_with_coords, summary_df, chart


    def create_coverage_plot(
        df: pl.DataFrame, width: int = 800, height: int = 300
    ) -> alt.Chart:
        """
        Create a coverage plot showing query coverage across targets

        Args:
            df: DataFrame with alignment data
            width: Chart width in pixels
            height: Chart height in pixels

        Returns:
            Altair chart showing coverage
        """
        coverage_data = df.group_by(["target_name", "query_name"]).agg(
            [
                pl.col("query_coverage").sum().alias("total_coverage"),
                pl.col("blast_identity").mean().alias("avg_identity"),
            ]
        )

        chart = (
            alt.Chart(coverage_data)
            .mark_bar()
            .encode(
                x=alt.X("target_name:N", title="Target Sequence"),
                y=alt.Y("total_coverage:Q", title="Query Coverage"),
                color=alt.Color(
                    "avg_identity:Q",
                    title="Avg Identity",
                    scale=alt.Scale(scheme="viridis"),
                ),
                tooltip=[
                    "target_name:N",
                    "query_name:N",
                    "total_coverage:Q",
                    "avg_identity:Q",
                ],
            )
            .properties(
                width=width, height=height, title="Query Coverage by Target"
            )
        )

        return chart
    return (load_and_visualize_paf,)


if __name__ == "__main__":
    app.run()
