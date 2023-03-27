import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt




def plot_islands(fst_file, namelist, overlapfile, plottitle):


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

    #getting overlap data
    with open(overlapfile, 'r') as comparison_overlap:
        overlaps = comparison_overlap.read().split('\n')


    #Create a new figure and set the size of the plot area within the figure
    fig, ax = plt.subplots(figsize=(60,5))
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

        #resetting the max
        if max(windowedfst) > max_fst:
            max_fst = max(windowedfst)

        #resetting min
        if min(windowedfst) < min_fst:
            min_fst = min(windowedfst)

    plt.axhline(y=pct1, linestyle='--')
    plt.axhline(y=pct2, linestyle='--')


    #Getting and plotting islands
    increment = 0.01

    for idx, comparison in enumerate(namelist):

        islands_data = islands[comparison].to_list()
        islands_data = [eval(x) if type(x) != float else x for x in islands_data]

        for i in islands_data:
            if isinstance(i, float):
                break
            else:
                start = i[0]
                stop = i[1]
                plt.hlines(y=max_fst+increment, xmin=float(start), xmax=float(stop), linewidth=15, color=plot_cols[idx])

        increment += 0.03

    for i in overlaps:
        if i == '':
            continue
        else:
            plt.fill_between(eval(i), min_fst, max_fst+increment, alpha=0.5, color='yellow')

    plt.title(plottitle)
    plt.xlabel("Position (bp)")
    plt.ylabel("Hudson Fst")
    plt.ticklabel_format(style='plain', axis='x', useOffset=False)
    for line in plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5)).get_lines():
        line.set_linewidth(4)
        line.set_alpha(1.0)

    plt.show()






plot_islands('win10kbp_ChrX_allcomparisons.csv', ['ZS.vs.RAL', 'ZS.vs.ZI'], 'win10kbp_ChrX_islands_ZS.vs.RAL_ZS.vs.ZI.txt', "Chrom X")
plot_islands('win10kbp_ChrX_allcomparisons.csv', ['ZS.vs.RAL', 'ZS.vs.FR', 'ZS.vs.ZI', 'ZS.vs.SAfr'], 'win10kbp_ChrX_islands_ZS.vs.RAL_ZS.vs.ZI_ZS.vs.FR_ZS.vs.SAfr.txt', "Chrom X")
plot_islands('win10kbp_ChrX_allcomparisons.csv', ['ZH.vs.RAL', 'ZH.vs.ZI'], 'win10kbp_ChrX_islands_ZH.vs.RAL_ZH.vs.ZI.txt', "Chrom X")
plot_islands('win10kbp_ChrX_allcomparisons.csv', ['ZW.vs.RAL', 'ZW.vs.ZI'], 'win10kbp_ChrX_islands_ZW.vs.RAL_ZW.vs.ZI.txt', "Chrom X")
plot_islands('win10kbp_ChrX_allcomparisons.csv', ['RAL.vs.FR', 'RAL.vs.ZI', 'RAL.vs.SAfr'], 'win10kbp_ChrX_islands_RAL.vs.FR_RAL.vs.ZI_RAL.vs.SAfr.txt', "Chrom X")



plot_islands('win10kbp_Chr2L_allcomparisons.csv', ['ZS.vs.RAL', 'ZS.vs.ZI'], 'win10kbp_Chr2L_islands_ZS.vs.RAL_ZS.vs.ZI.txt', "Chrom 2L")
plot_islands('win10kbp_Chr2L_allcomparisons.csv', ['ZS.vs.RAL', 'ZS.vs.FR', 'ZS.vs.ZI', 'ZS.vs.SAfr'], 'win10kbp_Chr2L_islands_ZS.vs.RAL_ZS.vs.ZI_ZS.vs.FR_ZS.vs.SAfr.txt', "Chrom 2L")
plot_islands('win10kbp_Chr2L_allcomparisons.csv', ['ZH.vs.RAL', 'ZH.vs.ZI'], 'win10kbp_Chr2L_islands_ZH.vs.RAL_ZH.vs.ZI.txt', "Chrom 2L")
plot_islands('win10kbp_Chr2L_allcomparisons.csv', ['ZW.vs.RAL', 'ZW.vs.ZI'], 'win10kbp_Chr2L_islands_ZW.vs.RAL_ZW.vs.ZI.txt', "Chrom 2L")
plot_islands('win10kbp_Chr2L_allcomparisons.csv', ['RAL.vs.FR', 'RAL.vs.ZI', 'RAL.vs.SAfr'], 'win10kbp_Chr2L_islands_RAL.vs.FR_RAL.vs.ZI_RAL.vs.SAfr.txt', "Chrom 2L")



plot_islands('win10kbp_Chr2R_allcomparisons.csv', ['ZS.vs.RAL', 'ZS.vs.ZI'], 'win10kbp_Chr2R_islands_ZS.vs.RAL_ZS.vs.ZI.txt', "Chrom 2R")
plot_islands('win10kbp_Chr2R_allcomparisons.csv', ['ZS.vs.RAL', 'ZS.vs.FR', 'ZS.vs.ZI', 'ZS.vs.SAfr'], 'win10kbp_Chr2R_islands_ZS.vs.RAL_ZS.vs.ZI_ZS.vs.FR_ZS.vs.SAfr.txt', "Chrom 2R")
plot_islands('win10kbp_Chr2R_allcomparisons.csv', ['ZH.vs.RAL', 'ZH.vs.ZI'], 'win10kbp_Chr2R_islands_ZH.vs.RAL_ZH.vs.ZI.txt', "Chrom 2R")
plot_islands('win10kbp_Chr2R_allcomparisons.csv', ['ZW.vs.RAL', 'ZW.vs.ZI'], 'win10kbp_Chr2R_islands_ZW.vs.RAL_ZW.vs.ZI.txt', "Chrom 2R")
plot_islands('win10kbp_Chr2R_allcomparisons.csv', ['RAL.vs.FR', 'RAL.vs.ZI', 'RAL.vs.SAfr'], 'win10kbp_Chr2R_islands_RAL.vs.FR_RAL.vs.ZI_RAL.vs.SAfr.txt', "Chrom 2R")



plot_islands('win10kbp_Chr3L_allcomparisons.csv', ['ZS.vs.RAL', 'ZS.vs.ZI'], 'win10kbp_Chr3L_islands_ZS.vs.RAL_ZS.vs.ZI.txt', "Chrom 3L")
plot_islands('win10kbp_Chr3L_allcomparisons.csv', ['ZS.vs.RAL', 'ZS.vs.FR', 'ZS.vs.ZI', 'ZS.vs.SAfr'], 'win10kbp_Chr3L_islands_ZS.vs.RAL_ZS.vs.ZI_ZS.vs.FR_ZS.vs.SAfr.txt', "Chrom 3L")
plot_islands('win10kbp_Chr3L_allcomparisons.csv', ['ZH.vs.RAL', 'ZH.vs.ZI'], 'win10kbp_Chr3L_islands_ZH.vs.RAL_ZH.vs.ZI.txt', "Chrom 3L")
plot_islands('win10kbp_Chr3L_allcomparisons.csv', ['ZW.vs.RAL', 'ZW.vs.ZI'], 'win10kbp_Chr3L_islands_ZW.vs.RAL_ZW.vs.ZI.txt', "Chrom 3L")
plot_islands('win10kbp_Chr3L_allcomparisons.csv', ['RAL.vs.FR', 'RAL.vs.ZI', 'RAL.vs.SAfr'], 'win10kbp_Chr3L_islands_RAL.vs.FR_RAL.vs.ZI_RAL.vs.SAfr.txt', "Chrom 3L")



plot_islands('win10kbp_Chr3R_allcomparisons.csv', ['ZS.vs.RAL', 'ZS.vs.ZI'], 'win10kbp_Chr3R_islands_ZS.vs.RAL_ZS.vs.ZI.txt', "Chrom 3R")
plot_islands('win10kbp_Chr3R_allcomparisons.csv', ['ZS.vs.RAL', 'ZS.vs.FR', 'ZS.vs.ZI', 'ZS.vs.SAfr'], 'win10kbp_Chr3R_islands_ZS.vs.RAL_ZS.vs.ZI_ZS.vs.FR_ZS.vs.SAfr.txt', "Chrom 3R")
plot_islands('win10kbp_Chr3R_allcomparisons.csv', ['ZH.vs.RAL', 'ZH.vs.ZI'], 'win10kbp_Chr3R_islands_ZH.vs.RAL_ZH.vs.ZI.txt', "Chrom 3R")
plot_islands('win10kbp_Chr3R_allcomparisons.csv', ['ZW.vs.RAL', 'ZW.vs.ZI'], 'win10kbp_Chr3R_islands_ZW.vs.RAL_ZW.vs.ZI.txt', "Chrom 3R")
plot_islands('win10kbp_Chr3R_allcomparisons.csv', ['RAL.vs.FR', 'RAL.vs.ZI', 'RAL.vs.SAfr'], 'win10kbp_Chr3R_islands_RAL.vs.FR_RAL.vs.ZI_RAL.vs.SAfr.txt', "Chrom 3R")