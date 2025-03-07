import pandas as pd
import numpy as np
import os


def range_len(ls):
    """
    Function takes a list with two elements in it representing a range and returns the difference between the two.
    Used in find_overlaps() to select the smallest range.
    """

    length = ls[1] - ls[0] 

    return length


def get_island_list(df, col):
    """
    Takes in the islands df and the column for interest and returns the column as a list of lists while removing nan values
    Returns islands as a list of 2-item lists
    """

    islands = df[col].to_list()
    islands = [eval(x) for x in islands if type(x) != float]

    return islands



def islands_above_value(island_csv_file, cutoff, chrom):

    print(island_csv_file)

    data = pd.read_csv(island_csv_file)
    cols = list(data.columns)

    
    for comp in cols:

        print(comp.replace("Win_FST_", "") + ":")

        comp_islands = get_island_list(data, comp)

        for island in comp_islands:
            length = range_len(island)

            if length > cutoff:
                island_str = str(island)
                island_str = island_str.replace(", ", "..")
                island_str = island_str.replace("[", chrom + ":")
                island_str = island_str.replace("]", "")

                print(island_str)

            else:
                continue
        
    print("\n")





islands_above_value("windowed_mean_5kbp_FST_averages_ChromX_islands.csv", 87850, "X")
print("\n")
islands_above_value("windowed_mean_5kbp_FST_averages_Chrom2L_islands.csv", 90000, "2L")
islands_above_value("windowed_mean_5kbp_FST_averages_Chrom2R_islands.csv", 90000, "2R")
islands_above_value("windowed_mean_5kbp_FST_averages_Chrom3L_islands.csv", 90000, "3L")
islands_above_value("windowed_mean_5kbp_FST_averages_Chrom3R_islands.csv", 90000, "3R")
print("\n")

islands_above_value("windowed_max_5kbp_FST_averages_ChromX_islands.csv", 83300, "X")
print("\n")
islands_above_value("windowed_max_5kbp_FST_averages_Chrom2L_islands.csv", 115000, "2L")
islands_above_value("windowed_max_5kbp_FST_averages_Chrom2R_islands.csv", 115000, "2R")
islands_above_value("windowed_max_5kbp_FST_averages_Chrom3L_islands.csv", 115000, "3L")
islands_above_value("windowed_max_5kbp_FST_averages_Chrom3R_islands.csv", 115000, "3R")
print("\n")