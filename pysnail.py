import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import SeqUtils
import os
import sys
import argparse
import numpy as np
from datetime import datetime

def read_fasta(filepath):
    data = []
    total_size = 0

    with open(filepath) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            GC = SeqUtils.GC(record.seq)
            length = len(record.seq)
            data.append((length, GC))
            total_size += length

    df = pd.DataFrame(data, columns = ["length", "GC"])

    df = df.sort_values("length", ascending = False)
    df = df.reset_index(drop = True)

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

def make_plot(df, total_size, largest, n50, n90, plotfile):
    offset = 0.5*np.pi
    theta = []
    theta_bg = []
    r = []
    gc_r = []
    r_bg = []
    r_gc_bg = []
    pos = 1

    bar_offset = offset + 2*np.pi*(largest/total_size)/2

    for i in range(1, total_size, 100000):
        theta_bg.append(offset + 2*np.pi/(total_size/i))
        r_bg.append(100)
        r_gc_bg.append(120)

    for row in df.iterrows():
        for i in range(pos, int(row[1]["cumlength"]), 100000):
            theta.append(offset + 2*np.pi/(total_size/i))
            r.append(100*((largest - row[1]["length"])/largest))
            gc_r.append(100 + row[1]["GC"]*0.2)
        pos += int(row[1]["length"])

    plt.figure(figsize = (10, 10))

    ax = plt.subplot(111, polar = True)

    plt.axis('off')

    ax.fill(theta_bg, r_gc_bg, color = "#a6cee3")
    ax.fill(theta, gc_r)

    ax.fill(theta_bg, r_bg, color = "#bbbbbb") # grey background of scaffolds
    ax.fill(theta, r, color = "#ffffff") # negative space of scaffolds

    scaffolds_bar = ax.bar(0, 0, color = "#bbbbbb")

    # plot the largest scaffold
    largest_bar = ax.bar(
        x=bar_offset, 
        height=df["length_rescale"][0], 
        width=df["width"][0], 
        bottom=100 - df["length_rescale"][0],
        color = "#e31a1c")

    # plot the n50 bar
    n50_bar = ax.bar(
        x = bar_offset + 2*np.pi*0.25 - max(df["width"])/2,
        height = 100*(n50/largest),
        width = 2*np.pi*0.5,
        bottom=100*(largest - n50)/largest,
        color="#ff7f00")

    # plot the n90 bar
    n90_bar = ax.bar(
        x = bar_offset + 2*np.pi*0.45 - max(df["width"])/2,
        height = 100*(n90/largest),
        width = 2*np.pi*0.9,
        bottom=100*(largest-n90)/largest,
        color="#fdbf6f")

    # various outlines added to plot
    ax.plot(theta, [100] * len(theta), color = "black")
    ax.plot(theta, [120] * len(theta), color = "black")
    ax.plot([offset, offset], [1, 100], color = "black")

    for i in np.arange(0, 360, 7.2):
        ax.plot([offset + np.radians(i), offset + np.radians(i)], [100, 102], color = "black")
        ax.plot([offset + np.radians(i), offset + np.radians(i)], [118, 120], color = "black")

    # apply percentage labels to GC axis
    for i in range(0, 360, 36):
        ax.text(offset + np.radians(i), 112, f"{int(100*(i/360))}%", rotation = i, ha = "center", va = "center", fontsize = 14)

    plt.title("Assembly statistics", fontsize = 20)

    plt.legend([scaffolds_bar, largest_bar, n50_bar, n90_bar], [
        f"Scaffold length (total {human_readable(total_size)})",
        f"Longest scaffold ({human_readable(largest)})", 
        f"N50 length ({human_readable(n50)})", 
        f"N90 length ({human_readable(n90)})"],
        loc=2,
        bbox_to_anchor=(0, 1.05),
        frameon = False)

    plt.savefig(plotfile)

def human_readable(base_length):
    i = 0
    units = [' bp', ' kB', ' MB', ' GB', ' TB']
    while base_length >= 1000:
        base_length = base_length / 1000
        i += 1
    return(f"{round(base_length, 1)}{units[i]}")

def main(arguments):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('infile', help="Input file", type=argparse.FileType('r'))
    parser.add_argument('-o', '--outfile', help="Output plot file",
                        default="plot.png", type=argparse.FileType('w'))

    args = parser.parse_args(arguments)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Hello! Loading fasta...")
    df, assembly_size, largest_scaffold = read_fasta(args.infile.name)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Total assembly size: {human_readable(assembly_size)}")
    
    n50, n90 = assembly_stats(df, assembly_size)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] n50: {human_readable(n50)} n90: {human_readable(n90)}")

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Making plot...")

    make_plot(df, assembly_size, largest_scaffold, n50, n90, args.outfile.name)

    print(f"[{datetime.now().strftime('%H:%M:%S')}] Plot finished! Bye!")

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))