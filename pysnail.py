import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import SeqUtils
import os
import sys
import argparse
import numpy as np
from datetime import datetime
import json


def read_fasta(filepath):
    data = []
    total_size = 0

    with open(filepath) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequence = record.seq.upper()
            GC = SeqUtils.GC(sequence)
            length = len(record.seq)
            N = 100*(sequence.count("N")/length)
            data.append((length, GC, N))
            total_size += length

    df = pd.DataFrame(data, columns=["length", "GC", "N"])

    df = df.sort_values("length", ascending=False)
    df = df.reset_index(drop=True)

    gaps = []
    for i in range(len(df["length"])):
        if i == 0:
            gap = 0
        else:
            gap = (df["length"][i] + df["length"][i-1]) / 2
        gaps.append(gap)

    df["gap"] = gaps

    largest = max(df["length"])
    df["length_rescale"] = 100 * (df["length"] / largest)
    df["width"] = 2*np.pi / (total_size/df["length"])
    df["angle"] = 2*np.pi / (total_size/df["gap"])
    df["angle"] = df["angle"].cumsum()
    df["cumlength"] = df["length"].cumsum()

    return df, total_size, largest


def assembly_stats(df, assembly_size):
    n50_sum = 0
    n90_sum = 0
    for length in df["length"]:
        if n50_sum < assembly_size*0.5:
            n50 = length
            n50_sum += n50
        if n90_sum < assembly_size*0.9:
            n90 = length
            n90_sum += n90

    return n50, n90


def make_plot(args, df, total_size, largest, n50, n90):
    offset = 0.5*np.pi # angle offset to make plot vertical
    theta = [] # angle for scaffold bins
    theta_bg = [] # angle for background bins (slightly different to theta, need to fix and merge)
    r = [] # size of scaffold for height, used to create negative space
    n_r = [] # size of N for scaffold
    gc_r = [] # size of gc for scaffold
    r_bg = [] # for background of scaffolds
    r_gc_bg = [] # gc background (ends up being AT)
    optional_offset = 135

    stepsize = total_size * 0.001 # to make sure plot doesn't degrade at low total sizes

    bar_offset = offset + 2*np.pi*(largest/total_size)/2

    for i in range(1, total_size, int(stepsize)):
        theta_bg.append(offset + 2*np.pi/(total_size/i))
        r_bg.append(100)
        r_gc_bg.append(120)

    pos = 1
    for row in df.iterrows():
        for i in range(pos, int(row[1]["cumlength"]), int(stepsize)):
            theta.append(offset + 2*np.pi/(total_size/i))
            r.append(100*((largest - row[1]["length"])/largest))
            n_r.append(100 + row[1]["N"]*0.2)
            gc_r.append(100 + row[1]["GC"]*0.2 + row[1]["N"]*0.2)
        pos += int(row[1]["length"])

    plt.figure(figsize=(10, 10))

    ax = plt.subplot(111, polar=True)

    plt.axis('off')

    ax.fill(theta_bg, r_gc_bg, color="#a6cee3") # gc background col
    ax.fill(theta, gc_r, color = "#1f78b4") # gc col
    ax.fill(theta, n_r, color = "#dddddd") # N col

    ax.fill(theta_bg, r_bg, color="#bbbbbb")  # grey background of scaffolds
    ax.fill(theta, r, color="#ffffff")  # negative space of scaffolds

    gc_bar = ax.bar(0, 0, color = "#1f78b4")
    at_bar = ax.bar(0, 0, color = "#a6cee3")
    scaffolds_bar = ax.bar(0, 0, color="#bbbbbb")

    # plot the largest scaffold
    largest_bar = ax.bar(
        x=bar_offset,
        height=df["length_rescale"][0],
        width=df["width"][0],
        bottom=100 - df["length_rescale"][0],
        color="#e31a1c")

    # plot the n50 bar
    n50_bar = ax.bar(
        x=bar_offset + 2*np.pi*0.25 - max(df["width"])/2,
        height=100*(n50/largest),
        width=2*np.pi*0.5,
        bottom=100*(largest - n50)/largest,
        color="#ff7f00")

    # plot the n90 bar
    n90_bar = ax.bar(
        x=bar_offset + 2*np.pi*0.45 - max(df["width"])/2,
        height=100*(n90/largest),
        width=2*np.pi*0.9,
        bottom=100*(largest-n90)/largest,
        color="#fdbf6f")

    # various outlines added to plot
    ax.plot(theta, [100] * len(theta), color="black")
    ax.plot(theta, [120] * len(theta), color="black")
    ax.plot([offset, offset], [1, 100], color="black")

    ax.plot(0.01, 140, color = None) # ghost to maintain plot size

    for i in np.arange(0, 360, 7.2):
        ax.plot([offset + np.radians(i), offset + np.radians(i)],
                [100, 102], color="black")
        ax.plot([offset + np.radians(i), offset + np.radians(i)],
                [118, 120], color="black")

    if args.busco is not None:
        
        busco_data = parse_busco(args.busco)

        lineage = busco_data["lineage_dataset"]["name"]
        busco_complete = busco_data["results"]["Complete"]
        busco_duplicate = busco_data["results"]["Multi copy"]
        busco_fragment = busco_data["results"]["Fragmented"]

        ax.bar(
            x = 1.6 * np.pi,
            align = "edge",
            height = 15,
            width = 0.4 * np.pi,
            bottom = optional_offset,
            color = "#dddddd")

        complete_bar = ax.bar(x = 1.6 * np.pi,
        align = "edge",
        height = 15,
        width = busco_complete/100 * 0.4 * np.pi,
        bottom = optional_offset,
        color = "#33a02c")

        duplicated_bar = ax.bar(x = 1.6 * np.pi,
        align = "edge",
        height = 15,
        width = busco_duplicate/100 * 0.4 * np.pi,
        bottom = optional_offset,
        color = "#236c1e")

        fragmented_bar = ax.bar(x = 1.6 * np.pi + (busco_complete/100 * 0.4 * np.pi), # offset by size of complete
        align = "edge",
        height = 15,
        width = busco_fragment/100 * 0.4 * np.pi,
        bottom = optional_offset,
        color = "#b2df8a")

        ax.bar(
            x = 1.6 * np.pi,
            align = "edge",
            height = 15,
            width = 0.4 * np.pi,
            bottom = optional_offset,
            fill = False,
            edgecolor = "black")

        busco_legend = plt.legend([complete_bar, duplicated_bar, fragmented_bar], [
        f"Complete ({busco_complete}%)",
        f"Duplicated ({busco_duplicate}%)",
        f"Fragmented ({busco_fragment}%)"],
        title = f"BUSCO {lineage}",
        loc = 4,
        frameon = False)

    if args.kmer is not None:
        # kmer completeness
        ax.bar(
            x = 0.01*np.pi,
            align = "edge",
            height = 15,
            width = 0.2 * np.pi,
            bottom = optional_offset,
            color = "#dddddd")

        kmer_bar = ax.bar(
            x = 0.01*np.pi,
            align = "edge",
            height = 15,
            width = args.kmer * 0.2 * np.pi,
            bottom = optional_offset,
            color = "#985f99")
        
        ax.bar(
            x = 0.01*np.pi,
            align = "edge",
            height = 15,
            width = 0.2 * np.pi,
            bottom = optional_offset,
            fill = False,
            edgecolor = "black")

        kmer_legend = plt.legend([kmer_bar], [f" K* ({args.kmer})"],
        title = "Completeness",
        loc = 1,
        frameon = False)

    # apply percentage labels to GC axis
    for i in range(0, 360, 36):
        ax.text(offset + np.radians(i), 112,
                f"{int(100*(i/360))}%", rotation=i, ha="center", va="center", fontsize=12)
        ax.text(offset + np.radians(i), 125,
                f"{human_readable(total_size*(i/360))}", rotation = i, ha = "center", va = "center", fontsize = 10)

    plt.title(args.title, fontsize=20)

    stats_legend = plt.legend([scaffolds_bar, largest_bar, n50_bar, n90_bar], [
        f"Scaffold length (total {human_readable(total_size)})",
        f"Longest scaffold ({human_readable(largest)})",
        f"N50 length ({human_readable(n50)})",
        f"N90 length ({human_readable(n90)})"],
        loc=2,
        frameon=False)

    ax = plt.gca().add_artist(stats_legend)

    bases_legend = plt.legend([gc_bar, at_bar], [
        f"GC (45.3%)",
        f"AT (54.7%)"],
        title = "Base composition",
        loc = 3,
        frameon = False)

    ax = plt.gca().add_artist(bases_legend)

    if args.kmer is not None:
        ax = plt.gca().add_artist(kmer_legend)

    if args.busco is not None:
        ax = plt.gca().add_artist(busco_legend)

    plt.savefig(args.outfile)


def human_readable(base_length):
    i = 0
    units = [' bp', ' kB', ' MB', ' GB', ' TB']
    while base_length >= 1000:
        base_length = base_length / 1000
        i += 1
    return (f"{round(base_length, 1)}{units[i]}")

def parse_busco(busco_json):
    infile = open(busco_json, "r")
    data = json.loads(infile.read())
    return data

def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return x


def main(arguments):

    parser = argparse.ArgumentParser(
        description = "A simple script to generate snail plots from fasta files",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('infile', help="Assembly fasta file for plot")
    parser.add_argument('-o', '--outfile',
                        help="Output plot file", default="plot.png")
    parser.add_argument('-t', '--title', help="Title of plot",
                        default="Assembly statistics")
    parser.add_argument("-b", "--busco", help = "BUSCO json file for score plotting")
    parser.add_argument("-k", "--kmer", type = restricted_float, help = "K* completeness value")

    args = parser.parse_args(arguments)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Hello! Loading fasta...")
    df, assembly_size, largest_scaffold = read_fasta(args.infile)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Total assembly size: {human_readable(assembly_size)}")

    n50, n90 = assembly_stats(df, assembly_size)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] n50: {human_readable(n50)} n90: {human_readable(n90)}")

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Making plot...")

    make_plot(args, df, assembly_size, largest_scaffold, n50, n90)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Plot finished! Bye!")


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
