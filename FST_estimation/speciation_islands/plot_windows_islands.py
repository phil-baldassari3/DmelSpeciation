import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt




def plot_islands(fst_file, namelist, plottitle):


    #opening fst data
    df = pd.read_csv(fst_file)

    #getting filenames for corresponding island and percentile files
    islandsname = fst_file.replace(".csv", "") + "_islands.csv"
    percentilename = fst_file.replace(".csv", "") + "_percentiles.txt"

    #opening island data
    islands = pd.read_csv(islandsname)

    #getting percentiles
    with open(percentilename, 'r') as pfile:
        data = pfile.read()
        lines = data.split('\n')

    pct1 = float(lines[0].split(' ')[2])
    pct2 = float(lines[1].split(' ')[2])

    #getting genomic position data
    pos = df['window_start'].to_list()

    """
    #getting overlap data
    with open(overlapfile, 'r') as comparison_overlap:
        overlaps = comparison_overlap.read().split('\n')
    """

    #Create a new figure and set the size of the plot area within the figure
    fig, ax = plt.subplots(figsize=(24,5))
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

    #Colors, will be used for later
    plot_cols = ['blue', 'orange', 'green', 'red', 'purple','brown', 'pink', 'gray', 'olive', 'cyan']

    #getting min and max so that islands dont overlap with fst (this is contnuously reset in the next block of code)
    min_fst = 0
    max_fst = 0

    #plotting fst data
    for idx, comparison in enumerate(namelist):

        windowedfst = df[comparison].to_list()
        plt.plot(pos, windowedfst, linewidth=0.5, alpha=0.5, color=plot_cols[idx], label=comparison)
        #plt.plot(pos, windowedfst, linewidth=1, alpha=0.8, color=plot_cols[idx], label=comparison)

        #resetting the max
        if max(windowedfst) > max_fst:
            max_fst = max(windowedfst)

        #resetting min
        if min(windowedfst) < min_fst:
            min_fst = min(windowedfst)

    plt.axhline(y=pct1, linestyle='--')
    plt.axhline(y=pct2, linestyle='--')


    #Getting and plotting islands
    increment = 0.007

    for idx, comparison in enumerate(namelist):

        islands_data = islands[comparison].to_list()
        islands_data = [eval(x) if type(x) != float else x for x in islands_data]

        for i in islands_data:
            if isinstance(i, float):
                break
            else:
                start = i[0]
                stop = i[1]
                plt.hlines(y=max_fst+increment, xmin=float(start), xmax=float(stop), linewidth=7, color=plot_cols[idx])

        increment += 0.03

    """
    for i in overlaps:
        if i == '':
            continue
        else:
            plt.fill_between(eval(i), min_fst, max_fst+increment, alpha=0.5, color='yellow')
    """

    plt.title(plottitle)
    plt.xlabel("Position (bp)")
    plt.ylabel("Hudson Fst")
    plt.ticklabel_format(style='plain', axis='x', useOffset=False)
    for line in plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5)).get_lines():
        line.set_linewidth(4)
        line.set_alpha(1.0)

    #plt.savefig(plottitle + ".png")
    plt.show()






plot_islands("windowed_mean_5kbp_FST_averages_ChromX.csv", ["Win_FST_Cos", "Win_FST_ZR", "Win_FST_ZH", "Win_FST_ZS"], "Chrom X Fst_fullwin")
plot_islands("windowed_mean_5kbp_FST_averages_Chrom2L.csv", ["Win_FST_Cos", "Win_FST_ZR", "Win_FST_ZH", "Win_FST_ZS"], "Chrom 2L Fst_fullwin")
plot_islands("windowed_mean_5kbp_FST_averages_Chrom2R.csv", ["Win_FST_Cos", "Win_FST_ZR", "Win_FST_ZH", "Win_FST_ZS"], "Chrom 2R Fst_fullwin")
plot_islands("windowed_mean_5kbp_FST_averages_Chrom3L.csv", ["Win_FST_Cos", "Win_FST_ZR", "Win_FST_ZH", "Win_FST_ZS"], "Chrom 3L Fst_fullwin")
plot_islands("windowed_mean_5kbp_FST_averages_Chrom3R.csv", ["Win_FST_Cos", "Win_FST_ZR", "Win_FST_ZH", "Win_FST_ZS"], "Chrom 3R Fst_fullwin")

plot_islands("windowed_max_5kbp_FST_averages_ChromX.csv", ["Win_FST_Cos", "Win_FST_ZR", "Win_FST_ZH", "Win_FST_ZS"], "Chrom X Fst_maxsnp")
plot_islands("windowed_max_5kbp_FST_averages_Chrom2L.csv", ["Win_FST_Cos", "Win_FST_ZR", "Win_FST_ZH", "Win_FST_ZS"], "Chrom 2L Fst_maxsnp")
plot_islands("windowed_max_5kbp_FST_averages_Chrom2R.csv", ["Win_FST_Cos", "Win_FST_ZR", "Win_FST_ZH", "Win_FST_ZS"], "Chrom 2R Fst_maxsnp")
plot_islands("windowed_max_5kbp_FST_averages_Chrom3L.csv", ["Win_FST_Cos", "Win_FST_ZR", "Win_FST_ZH", "Win_FST_ZS"], "Chrom 3L Fst_maxsnp")
plot_islands("windowed_max_5kbp_FST_averages_Chrom3R.csv", ["Win_FST_Cos", "Win_FST_ZR", "Win_FST_ZH", "Win_FST_ZS"], "Chrom 3R Fst_maxsnp")


