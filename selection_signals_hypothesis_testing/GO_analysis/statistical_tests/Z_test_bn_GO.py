import pandas as pd
import numpy as np
import random
import scipy.stats
from statistics import mean, stdev
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.ticker import MaxNLocator
import seaborn as sns
import warnings
warnings.filterwarnings("ignore", message="use_inf_as_na option is deprecated")


def permutation_Ztest(test_set, null_set, permutations, direction, popname, GOterm):
    """
    Function performs a Z-test to compare the mean genomic statistic of the GO category of interest
    to the mean genomic statistic of Non-Go category genes. Permutating is done to downsample the
    the null set so that the sample size between the two sets are the same for the comparison

    Arguments:
    test_set (list): list of genomic statistic values for the GO category of interest
    null_set (list): list of genomic statistic values for the Non-GO category of interest
    permutations (int): number of sampling permutations for the null set
    direction (str): "left", "right", or "two-tailed"

    Returns test_avg, sampled_null_list_of_avgs, zscore, pvalue
    """

    #finding sample size and test average
    sample_size = len(test_set)
    test_avg = mean(test_set)

    #permutating
    sampled_null_list_of_avgs = []

    with tqdm(range(permutations), desc=f"Permutating for {popname} & {GOterm}") as pbar:
        for i in range(permutations):
            null_sample = random.sample(null_set, sample_size)
            null_sample_avg = mean(null_sample)
            sampled_null_list_of_avgs.append(null_sample_avg)

            pbar.update()

    #computing Z-score
    zscore = (test_avg - mean(sampled_null_list_of_avgs)) / stdev(sampled_null_list_of_avgs)

    #compute p-value
    if direction == "right":
        pvalue = scipy.stats.norm.sf(zscore)
    elif direction == "left":
        pvalue = scipy.stats.norm.cdf(zscore)
    elif direction == "two-tailed":
        if zscore < 0:
            pvalue = 2 * (scipy.stats.norm.cdf(zscore))
        else:
            pvalue = 2 * (scipy.stats.norm.sf(zscore))
    else:
        raise(SyntaxError, "Wrong argument given for `direction`")
    

    return test_avg, sampled_null_list_of_avgs, zscore, pvalue




def parse_data_run_permutation_tests(datacsv, popnames, GOterms, permutations, direction):
    """
    Function parses data to run the Ztest for each population and GO term.

    Arguments:
    data (str): filename for genomic statistic data csv
    popnames (list): list of pop names to run tests on, must match csv column names
    GOterms (list): list of GO terms to run tests on, must match csv column names
    permutations (int): number of sampling permutations for the null set
    direction (str): "left", "right", or "two-tailed"

    Returns df_of_GOavgs, df_of_null_sample_lists, df_of_zcores, df_of_pvalues
    """

    #loading data
    data = pd.read_csv(datacsv)

    #setting empty dfs
    df_of_GOavgs = pd.DataFrame()
    df_of_null_sample_lists = pd.DataFrame()
    df_of_zcores = pd.DataFrame()
    df_of_pvalues = pd.DataFrame()

    df_of_GOavgs["Population"] = popnames
    df_of_null_sample_lists["Population"] = popnames
    df_of_zcores["Population"] = popnames
    df_of_pvalues["Population"] = popnames


    #looping through pops and GOs
    for go in GOterms:

        #setting empty lists that will become df column
        col_list_GOavgs = []
        col_list_nulls = []
        col_list_zscores = []
        col_list_pvalues = []

        for pop in popnames:
            
            #filtering and grabbing data
            GOdf = data[data[go] == go.replace("_", " ")]
            NonGOdf = data[data[go] != go.replace("_", " ")]

            GOset = GOdf[pop].to_list()
            NonGOset = NonGOdf[pop].to_list()

            #run Z-test
            GOavg, NonGOavgs, z_score, p_value = permutation_Ztest(GOset, NonGOset, permutations, direction, pop, go)

            #appending to col
            col_list_GOavgs.append(GOavg)
            col_list_nulls.append(NonGOavgs)
            col_list_zscores.append(z_score)
            col_list_pvalues.append(p_value)

        #adding column to dfs
        df_of_GOavgs[go] = col_list_GOavgs
        df_of_null_sample_lists[go] = col_list_nulls
        df_of_zcores[go] = col_list_zscores
        df_of_pvalues[go] = col_list_pvalues


    return df_of_GOavgs, df_of_null_sample_lists, df_of_zcores, df_of_pvalues



def Zscore_heatmap(table, title):
    """
    function plots a heatmap from the Z-score table generated from parse_data_run_permutation_tests.
    """

    #setting row labels to index
    table = table.set_index(table.columns[0])
    table.columns = table.columns.str.replace("_", " ", regex=False)

    #determine absolute max
    abs_max = max(abs(table.values.min()), abs(table.values.max()))

    #setting color map
    cmap = plt.cm.seismic

    #plot heat
    plt.figure(figsize=(10, 10))
    heat = plt.imshow(table.values, cmap=cmap, vmin=-abs_max, vmax=abs_max)

    #adding scale bar
    cbar = plt.colorbar(heat, shrink=0.8, pad=0.02)
    cbar.set_label("Z-score", rotation=270, labelpad=15)

    #adding gridlines
    num_rows, num_cols = table.shape
    for r in range(1, num_rows):
        plt.axhline(r - 0.5, color="white", linewidth=1)
    for c in range(1, num_cols):
        plt.axvline(c - 0.5, color="white", linewidth=1)

    #adding labels
    plt.xticks(range(table.shape[1]), table.columns, rotation=45, ha="right")
    plt.yticks(range(table.shape[0]), table.index)

    #formatting plots
    plt.grid(False)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.title(title)
    plt.subplots_adjust(left=0.2, bottom=0.2)
    #plt.show()

    plt.savefig(f"{title}_zcores.png", dpi=300, bbox_inches="tight")
    plt.close()


def Pvalue_heatmap(table, title):
    """
    function plots a heatmap from the p-value table generated from parse_data_run_permutation_tests.
    """

    #setting row labels to index
    table = table.set_index(table.columns[0])
    table.columns = table.columns.str.replace("_", " ", regex=False)

    #setting color map
    bounds = [0, 0.001, 0.01, 0.05, 1.0]
    colors = ["darkred", "red", "pink", "lightgray"]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(bounds, cmap.N)

    #plot heat
    plt.figure(figsize=(10, 10))
    heat = plt.imshow(table.values, cmap=cmap, norm=norm)

    #adding scale bar
    cbar = plt.colorbar(heat, shrink=0.8, pad=0.02, ticks=[0.0005, 0.005, 0.03, 0.75])
    cbar.ax.set_yticklabels(["p≤0.001", "0.001<p≤0.01", "0.01<p≤0.05", "p>0.05"])

    #adding gridlines
    num_rows, num_cols = table.shape
    for r in range(1, num_rows):
        plt.axhline(r - 0.5, color="white", linewidth=1)
    for c in range(1, num_cols):
        plt.axvline(c - 0.5, color="white", linewidth=1)

    #adding labels
    plt.xticks(range(table.shape[1]), table.columns, rotation=45, ha="right")
    plt.yticks(range(table.shape[0]), table.index)

    #formatting plots
    plt.grid(False)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.title(title)
    plt.subplots_adjust(left=0.2, bottom=0.2)
    #plt.show()

    plt.savefig(f"{title}_pvalues.png", dpi=300, bbox_inches="tight")
    plt.close()


def plot_nulls_and_avgs(avgs_table, nulls_table, statistic, title, color):
    """Function makes a faceted plot showing the null distibutions and the GO average for each Pop/GO combination"""

    #Melt the data into long format
    avg_melted = avgs_table.melt(id_vars=avgs_table.columns[0], var_name="GO", value_name="avg")
    null_melted = nulls_table.melt(id_vars=nulls_table.columns[0], var_name="GO", value_name="null_list")

    #Merge the long dfs
    merged = pd.merge(avg_melted, null_melted, on=["Population", "GO"])

    #Explode the null distibution lists
    merged_exploded = merged.explode("null_list")
    merged_exploded["null_list"] = pd.to_numeric(merged_exploded["null_list"], errors='coerce')

    #clean GO labels
    merged_exploded["GO"] = merged_exploded["GO"].str.replace("_", " ")

    # Set up the facet grid
    g = sns.FacetGrid(merged_exploded, row="Population", col="GO", margin_titles=True, sharex=False, sharey=False, height=2.5, aspect=1.5, despine=False)

    # Plot the histogram
    g.map_dataframe(sns.histplot, x="null_list", color="black", bins=50, stat="density")

    # Add vertical line for GO average
    def draw_vline(data, color, **kwargs):
        plt.axvline(x=data["avg"].iloc[0], color=color, linewidth=2.5)

    g.map_dataframe(draw_vline, color=color)

    #format labels
    g.set_titles(row_template="{row_name}", col_template="{col_name}", size=12, weight="bold")

    for ax in g.axes.flat:
        #ax.title.set_fontsize(12)
        #ax.title.set_fontweight("bold")
        ax.yaxis.set_ticks([])
        ax.set_yticklabels([])
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, prune='both'))


    #gloabl labels
    g.set_axis_labels("", "")

    g.figure.text(0.5, 0.04, statistic, ha="center", va="center", fontsize=12)
    g.figure.text(0.04, 0.5, "Density", ha="center", va="center", rotation="vertical", fontsize=12)
    g.figure.suptitle(title, fontsize=14)
    g.tight_layout(rect=[0.05, 0.05, 1, 0.95])


    #plt.show()

    plt.savefig(f"{title}_permutation.png", dpi=300, bbox_inches="tight")
    plt.close()









random.seed(761)

poplist = ["ZS", "ZH", "RAL", "FR", "ZI"]
popcomplist = ["ZS.vs.RAL", "ZS.vs.FR", "ZS.vs.ZI", "ZH.vs.RAL", "ZH.vs.FR", "ZH.vs.ZI", "RAL.vs.FR", "FR.vs.ZI", "ZI.vs.RAL"]
GOlist = ["nervous_system_development", "sensory_system_development", "mechanosensory_behavior", "aggressive_behavior", "learning_or_memory", "nervous_system_process"]



def main_func(csvfile, population_label_list, GOterm_list, num_permutations, significance_direction, genomic_stat, title, color):
    """Main Function to run the permuation test and output plots"""

    avgsDF, nullDF, zDF, pDF = parse_data_run_permutation_tests(csvfile, population_label_list, GOterm_list, num_permutations, significance_direction)

    #Zscore_heatmap(zDF, title)
    #Pvalue_heatmap(pDF, title)
    plot_nulls_and_avgs(avgsDF, nullDF, genomic_stat, title, color)

    pDF.to_csv(f"{title}_p-values.csv", index=False)









#main_func calls
main_func("GO_ann_FST_gene.csv", popcomplist, GOlist, 5000, "right", "FST", "FST (Gene-level)", "turquoise")
main_func("GO_ann_FST_win.csv", popcomplist, GOlist, 5000, "right", "FST", "FST (10kbp Windows)", "turquoise")
main_func("GO_ann_PI_gene.csv", poplist, GOlist, 5000, "left", "Pi", "Pi (Gene-level)", "orange")
main_func("GO_ann_PI_win.csv", poplist, GOlist, 5000, "left", "Pi", "Pi (10kbp Windows)", "orange")
main_func("GO_ann_THETA.csv", poplist, GOlist, 5000, "left", "Theta", "Watterson's Theta (10kbp Windows)", "crimson")
main_func("GO_ann_TjD.csv", poplist, GOlist, 5000, "left", "Tajima's D", "Tajima's D (10kbp Windows)", "purple")


