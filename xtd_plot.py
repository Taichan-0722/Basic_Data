import os
import sys
import re
import io
import csv
import glob
import subprocess
import shutil
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LogNorm

title = "Radius vs Norm"
xlabel = "Distance from the center (arcmin)"
ylabel = "Norm (/400$\pi$ arcmin$^2$)"
xlog = "off"
ylog = "off"
jitter_width = 0.1
point_format = "o"
point_formats = ["o", "^", "s", "v", "D"]
point_color = "c"
point_edge_color = "k"
point_edge_width = 0.5
errorbar_color = "grey"
errorbar_width = 1
errorbar_edge_length = 3
colorbar_min = 0
colorbar_max = 360
colorbar_type = "hsv"
colorbar_label = "Azimuth (deg) (counterclockwise starting from north)"

csv_data = input("\n Input the csv data name. > ").strip()

if not os.path.exists(csv_data):
    print(f"\n '{csv_data}' does'nt exist.\n")
    sys.exit(1)

jitter_or_not = input(" Create a graph with jitter? (yes|no) > ").strip()

colorbar_or_not = input(" Create the graph with the colorbar? (yes|no) > ").strip()

if colorbar_or_not == "yes":
    group_or_not = input(f" Group the points (up to {len(point_formats)})? (yes|no) > ").strip()

elif colorbar_or_not == "no":
    group_or_not = input(" Group the points (up to 20)? (yes|no) > ").strip()

if colorbar_or_not == "yes" or group_or_not == "yes":
    errorbar_color_or_not = input(" Match the error bar color to the point color? (yes|no) > ").strip()

else:
    errorbar_color_or_not = "no"

print()

df = pd.read_csv(csv_data)
norm = Normalize(vmin=colorbar_min, vmax=colorbar_max)
cmap = plt.cm.get_cmap(colorbar_type)

plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)

if xlog == "on":
    plt.xscale('log')

if ylog == "on":
    plt.yscale('log')

if colorbar_or_not == "yes" and group_or_not == "yes":

    groups = df["group"].unique()

    for i, grp in enumerate(groups):

        subset = df[df["group"] == grp]
        x = subset["x"]
        y = subset["y"]
        colorbar_values = subset["colorbar"]
        xerr = [subset["xerr_low"], subset["xerr_high"]]
        yerr = [subset["yerr_low"], subset["yerr_high"]]

        if jitter_or_not == "yes":
            x = x + np.random.uniform(-jitter_width, jitter_width, size=len(x))

        if errorbar_color_or_not == "yes":
            for xi, yi, xerr_low_i, xerr_high_i, yerr_low_i, yerr_high_i, colorbar_value_i in zip(x, y, xerr[0], xerr[1], yerr[0], yerr[1], colorbar_values):
                color_i = cmap(norm(colorbar_value_i))
                plt.errorbar([xi], [yi], xerr=[[xerr_low_i], [xerr_high_i]], yerr=[[yerr_low_i], [yerr_high_i]], fmt='none', ecolor=color_i, elinewidth=errorbar_width, capsize=errorbar_edge_length, zorder = 1)

        elif errorbar_color_or_not == "no":
            plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='none', ecolor=errorbar_color, elinewidth=errorbar_width, capsize=errorbar_edge_length, zorder = 1)

        sc = plt.scatter(x, y, marker=point_formats[i], linewidths=point_edge_width, edgecolors=point_edge_color, c=colorbar_values, cmap=colorbar_type, vmin=colorbar_min, vmax=colorbar_max, label=grp, zorder = 2)

    colorbar = plt.colorbar(sc, pad=0.02)
    colorbar.set_label(colorbar_label)
    colorbar.ax.tick_params(which="both", direction="out")
    plt.legend()
    plt.show()

elif colorbar_or_not == "yes" and group_or_not == "no":

    x = df["x"]
    y = df["y"]
    colorbar_values = df["colorbar"]
    xerr = [df["xerr_low"], df["xerr_high"]]
    yerr = [df["yerr_low"], df["yerr_high"]]

    if jitter_or_not == "yes":
        x = x + np.random.uniform(-jitter_width, jitter_width, size=len(x))

    if errorbar_color_or_not == "yes":
        for xi, yi, xerr_low_i, xerr_high_i, yerr_low_i, yerr_high_i, colorbar_value_i in zip(x, y, xerr[0], xerr[1], yerr[0], yerr[1], colorbar_values):
            color_i = cmap(norm(colorbar_value_i))
            plt.errorbar([xi], [yi], xerr=[[xerr_low_i], [xerr_high_i]], yerr=[[yerr_low_i], [yerr_high_i]], fmt='none', ecolor=color_i, elinewidth=errorbar_width, capsize=errorbar_edge_length, zorder = 1)

    elif errorbar_color_or_not == "no":
        plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='none', ecolor=errorbar_color, elinewidth=errorbar_width, capsize=errorbar_edge_length, zorder = 1)

    sc = plt.scatter(x, y, marker=point_format, linewidths=point_edge_width, edgecolors=point_edge_color, c=colorbar_values, cmap=colorbar_type, vmin=colorbar_min, vmax=colorbar_max, zorder = 2)
    colorbar = plt.colorbar(sc, pad=0.02)
    colorbar.set_label(colorbar_label)
    colorbar.ax.tick_params(which="both", direction="out")

    plt.show()

elif colorbar_or_not == "no" and group_or_not == "yes":

    groups = df["group"].unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(groups)))

    for i, grp in enumerate(groups):

        subset = df[df["group"] == grp]
        x = subset["x"]
        y = subset["y"]
        xerr = [subset["xerr_low"], subset["xerr_high"]]
        yerr = [subset["yerr_low"], subset["yerr_high"]]

        if jitter_or_not == "yes":
            x = x + np.random.uniform(-jitter_width, jitter_width, size=len(x))

        if errorbar_color_or_not == "yes":
            plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=point_format, mec=point_edge_color, mew=point_edge_width, ecolor=colors[i], elinewidth=errorbar_width, capsize=errorbar_edge_length, color=colors[i], label=grp)

        elif errorbar_color_or_not == "no":
            plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=point_format, mec=point_edge_color, mew=point_edge_width, ecolor=errorbar_color, elinewidth=errorbar_width, capsize=errorbar_edge_length, color=colors[i], label=grp)

    plt.legend()
    plt.show()

elif colorbar_or_not == "no" and group_or_not == "no":

    x = df["x"]
    y = df["y"]
    xerr = [df["xerr_low"], df["xerr_high"]]
    yerr = [df["yerr_low"], df["yerr_high"]]

    if jitter_or_not == "yes":
        x = x + np.random.uniform(-jitter_width, jitter_width, size=len(x))

    plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=point_format, mec=point_edge_color, mew=point_edge_width, ecolor=errorbar_color, elinewidth=errorbar_width, capsize=errorbar_edge_length, color=point_color)

    plt.show()