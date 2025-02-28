#! /usr/bin/env python

import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
import argparse
import os

def setup_and_parse_args():
    parser = argparse.ArgumentParser(description="generate the report of downsample.")
    parser.add_argument("-s", "--samples", required=True, help="sample names")
    parser.add_argument("-t", "--saturation_template", required=True, help="html template for saturation")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to the output path")
    args = parser.parse_args()
    return args

def plot(df, sample):
    color_palette = ["#5963f5", "#e44c39", "#37c58d", "#9c59f5"]

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=df["Downsample Ratio"],
        y=df["Sequencing Saturation"],
        name="Sequencing Saturation",
        yaxis="y3",
        line=dict(color=color_palette[0]),
    ))

    fig.add_trace(go.Scatter(
        x=df["Downsample Ratio"],
        y=df["Duplication Ratio"],
        name="Duplication Ratio",
        yaxis="y3",
        line=dict(color=color_palette[1]),
    ))

    fig.add_trace(go.Scatter(
        x=df["Downsample Ratio"],
        y=df["UMI Counts"],
        name="UMI Counts",
        yaxis="y2",
        line=dict(color=color_palette[2]),
    ))

    fig.add_trace(go.Scatter(
            x=df["Downsample Ratio"],
            y=df["UMI Types"],
            name="UMI Types",
            yaxis="y",
            line=dict(color=color_palette[3]),
        ))

    # style all the traces
    fig.update_traces(
        hoverinfo="y",
        mode="lines",
    )

    # Update axes
    fig.update_layout(
        title=dict(text=f"{sample}",
                   x=0.5
        ),
        xaxis=dict(
            title=dict(
                text="Downsample Ratio",
            ),
            title_standoff=0,
        ),
        yaxis=dict(
            anchor="x",
            autorange=True,
            domain=[0.05, 0.5],
            range=[0, None],
            showline=True,
            linecolor=color_palette[3],
            linewidth=1.5,
            title=dict(
                text="UMI Types",
                font=dict(color=color_palette[3])
            ),
            tickfont=dict(color=color_palette[3]),
            title_standoff=0
        ),
        yaxis2=dict(
            anchor="x",
            autorange=True,
            domain=[0.05, 0.5],
            range=[0, None],
            showline=True,
            overlaying="y",
            side="right",
            linecolor=color_palette[2],
            linewidth=1.5,
            title=dict(
                text="UMI Counts",
                font=dict(color=color_palette[2])
            ),
            tickfont=dict(color=color_palette[2]),
            title_standoff=0
        ),
        yaxis3=dict(
            anchor="x",
            domain=[0.55, 1],
            range=[0, 100],
            showline=True,
            linecolor='black',
            linewidth=1.5,
            title=dict(
                text="Saturation/%",
            ),
            title_standoff=0,
            dtick=20
        ),
        height=500,
        width=600,
        template="plotly_white",
        legend=dict(
                orientation='h',
                x=0.5,
                xanchor='center',
                y=-0.1,
                traceorder='normal',
                itemwidth=30
            ),
        legend_tracegroupgap=5
    )
    return fig

if __name__ == "__main__":
    args = setup_and_parse_args()
    samples = args.samples.split()
    saturation_template = args.saturation_template
    output_dir = args.output_dir
    summary_dir = os.path.join(output_dir, "00_summary")
    
    fig_html_list = []
    for sample in samples:
        saturation_dir = os.path.join(output_dir, sample, "04_saturation")
        file = os.path.join(saturation_dir, f"{sample}_Downsample.tsv")
        df = pd.read_csv(file, sep="\t")
        fig = plot(df, sample)
        fig_html = fig.to_html(full_html=False, include_plotlyjs='cdn')
        fig.write_html(f"{saturation_dir}/{sample}_Downsample.html", auto_open=False, full_html=True, include_plotlyjs='cdn')
        pio.full_figure_for_development(fig, warn=False)
        fig.write_image(f"{saturation_dir}/{sample}_Downsample.pdf", format='pdf')
        fig_html_list.append(fig_html)

    figs_html = "\n".join(f'<div class="plot">{fig_html}</div>' for fig_html in fig_html_list)

    with open(saturation_template, "r", encoding="utf-8") as f:
        saturation_template = f.read()

    final_html = saturation_template.replace("{{ figs_html }}", figs_html)
    with open(f"{summary_dir}/Saturation.html", "w", encoding="utf-8") as f:
        f.write(final_html)