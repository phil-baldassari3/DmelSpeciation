import pandas as pd
import numpy as np
from statistics import mean
import matplotlib.pyplot as plt



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


def get_island_lengths(islandls):
    """
    Takes in the island list and computes the length of each island
    Returns list of island lengths
    """

    island_lens = []
    for i in islandls:
        l = range_len(i)
        island_lens.append(l)

    return island_lens


def get_pct_of_lengths(island_csv_file, pct):

    data = pd.read_csv(island_csv_file)
    cols = list(data.columns)

    #get list of all island lengths
    all_lens = []
    for comp in cols:
        comp_islands = get_island_list(data, comp)
        comp_is_lens = get_island_lengths(comp_islands)
        all_lens += comp_is_lens

    #percentile
    lenpct = np.percentile(all_lens, pct)

    #filtering for lengths above percentile
    pct_lens = [i for i in all_lens if i > lenpct]

    #count regions above the percentile
    numregions = len(pct_lens)

    return lenpct, numregions






def make_summ_stat_table(island_csv_file):
    """
    Takes in an speciation islands csv file and outputs a csv that computes the
    1) total number of islands per comparisons
    2) cummulative length of islands per comparison
    3) average length of islands per comparison
    """

    data = pd.read_csv(island_csv_file)
    cols = list(data.columns)

    comps = []
    counts = []
    cummlens = []
    avglens = []

    for comp in cols:
        comp_islands = get_island_list(data, comp)
        comp_is_lens = get_island_lengths(comp_islands)



        """
        #optional plot the histogram
        plt.hist(comp_is_lens, bins=200)
        plt.xlabel('length')
        plt.ylabel('Frequency')
        plt.title(comp)
        plt.show()
        """


   
        comp_count = len(comp_islands)
        comp_cummlen = sum(comp_is_lens)
        if comp_count == 0:
            comp_avglen = 0
        else:
            comp_avglen = mean(comp_is_lens)

        comps.append(comp)
        counts.append(comp_count)
        cummlens.append(comp_cummlen)
        avglens.append(comp_avglen)


    dictionary = {"Comparison":comps, "Num_of_islands":counts, "Cummulative_len_of_islands":cummlens, "Average_island_len":avglens}

    summdf = pd.DataFrame(dictionary)

    summdf.to_csv("summary_stats_" + island_csv_file, index=False)





"""
make_summ_stat_table("windowed_mean_5kbp_FST_averages_ChromX_islands.csv")
make_summ_stat_table("windowed_mean_5kbp_FST_averages_Autosomes_islands.csv")
make_summ_stat_table("windowed_max_5kbp_FST_averages_ChromX_islands.csv")
make_summ_stat_table("windowed_max_5kbp_FST_averages_Autosomes_islands.csv")
"""


print(get_pct_of_lengths("windowed_mean_5kbp_FST_averages_ChromX_islands.csv", 99))
print(get_pct_of_lengths("windowed_mean_5kbp_FST_averages_Autosomes_islands.csv", 99))
print(get_pct_of_lengths("windowed_max_5kbp_FST_averages_ChromX_islands.csv", 99))
print(get_pct_of_lengths("windowed_max_5kbp_FST_averages_Autosomes_islands.csv", 99))