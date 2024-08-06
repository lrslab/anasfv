#!/usr/bin/env python
# coding: utf-8
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO
from pycirclize import Circos
from matplotlib.lines import Line2D
import os
import pandas as pd


def recombination_plot(file_name,format='pdf'):
    df = pd.read_csv(f'{file_name}', sep='\t')
    genome_len = df.End.iloc[-1]

    typeI_features = []
    typeII_features = []
    unknown_features = []

    for i in range(df.shape[0]):
        if df['Conclusion'][i] == 'I':
            typeI_features.append(SeqFeature(FeatureLocation(int(df.Start[i]), int(df.End[i])), type="typeI"))
        elif df['Conclusion'][i] == 'II':
            typeII_features.append(SeqFeature(FeatureLocation(int(df.Start[i]), int(df.End[i])), type="typeII"))
        else:
            unknown_features.append(SeqFeature(FeatureLocation(int(df.Start[i]), int(df.End[i])), type="unknown"))
    circos = Circos(sectors={'ASFV': genome_len})
    sector = circos.get_sector('ASFV')
    major_ticks_interval = 8000
    minor_ticks_interval = 1000
    outer_track = sector.add_track((73, 75))
    outer_track.axis(fc="lightgrey")
    outer_track.xticks_by_interval(
        major_ticks_interval, label_formatter=lambda v: f"{int(v / 10 ** 3)} kb"
    )
    outer_track.xticks_by_interval(minor_ticks_interval, tick_length=1, show_label=False)

    typeI_track = sector.add_track((63, 70), r_pad_ratio=0.1)
    typeI_track.genomic_features(typeI_features, fc="blue")
    typeII_track = sector.add_track((63, 70), r_pad_ratio=0.1)
    typeII_track.genomic_features(typeII_features, fc="red")
    unknown_track = sector.add_track((63, 70), r_pad_ratio=0.1)
    unknown_track.genomic_features(unknown_features, fc="gray")

    fig = circos.plotfig()

    # Add legend
    handles = [
        Line2D([], [], color="blue", label="Genotype I ASFV origin", marker='s', ms=8, ls="None"),
        Line2D([], [], color="red", label="Genotype II ASFV origin", marker='s', ms=8, ls="None"),
        Line2D([], [], color="gray", label="Uncertain ASFV origin", marker='s', ms=8, ls="None"),

    ]
    _ = circos.ax.legend(handles=handles, bbox_to_anchor=(0.5, 0.475), loc="center", fontsize=8)
    if file_name.endswith('.tsv'):
        file_name=file_name[:-4]
    fig.savefig(f"{file_name}.{format}")


if __name__ =='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='tsv file from recombination_test.py')
    parser.add_argument("-f", "--format", help="output format (jpg, png or pdf)",default='pdf')
    args = parser.parse_args()
    input_path = args.input
    format = args.format

    recombination_plot(input_path, format)