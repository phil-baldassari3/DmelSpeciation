import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def load_csv(csvfile):
    """
    Loads in the csv of windowed Fst data and returns a dataframe and a list of certain columns that will be used in the other functions
    """

    #loading dataframe
    df = pd.read_csv(csvfile)

    #listing of FST columns
    fst_cols = list(df.columns)
    fst_cols.remove("Chrom")
    fst_cols.remove("window_start")
    fst_cols.remove("window_end")
    fst_cols.remove("FBgn")
    fst_cols.remove("gene_symbol")

    return df, fst_cols


def high_FST4comparison(df, col, percentile):
    """
    Function takes the FST dataframe, the column name of interest, and a percentile (e.g. 95th pct is 95)
    The function finds the percentile cutoff of the column of interest, filters the dataframe rows to be above
    that cutoff, and returns the filtered dataframe, a list of FST values from the column of interest, and the
    cutoff value
    """

    #removing windows where all FST is zero
    df = df.loc[(df.iloc[:, 3:9] != 0).any(axis=1)]
    
    #getting FST list
    fst_ls = df[col].to_list()

    #get pct cutoff
    pct = np.percentile(fst_ls, percentile)

    #filtering df
    filtered_df = df[df[col] > pct]
    filtered_df = filtered_df[["Chrom", "window_start", "window_end", col, "FBgn", "gene_symbol"]]
    filtered_df = filtered_df.sort_values(by=[col], ascending=False)

    return filtered_df, fst_ls, pct


def plotting_dist(fst_ls, pctcutoff, title):
    """
    Plots the FST distribution with the percentile cutoff
    """

     #number of bins, Freedman–Diaconis
    q25, q75 = np.percentile(fst_ls, [25, 75])
    bin_width = 2 * (q75 - q25) * len(fst_ls) ** (-1/3)
    bins = round((max(fst_ls) - min(fst_ls)) / bin_width)

    #plotting
    plt.figure(figsize=(10, 8))
    plt.hist(fst_ls, density=True, bins=bins)
    plt.axvline(x=pctcutoff, color='r', lw=1, linestyle='--')
    plt.ylabel('Density')
    plt.xlabel("FST")
    plt.title(title)

    plt.savefig(f"{title}_dist.png")
    plt.close()


def output_files(filtered_df, name):
    """
    Function outputs the filtered dataframe as a csv and outputs a text file with a list of genes from the dataframe
    """

    #saving df
    filtered_df.to_csv(name + ".csv", index=False)

    #getting list of genes
    filtered_df.dropna(subset=["FBgn"], inplace=True)
    tempgenels = filtered_df["FBgn"].to_list()
    genestr = "; ".join(tempgenels)
    gene_ls = genestr.split("; ")

    genes = list(set(gene_ls))

    #saving gene list
    with open(name + "_genes.txt", 'w') as genes_file:
        genes_file.write("\n".join(genes))





def main_func(csv, pct):

    print(csv)

    output_name = "95thpct_" + csv.replace(".csv", "")

    data, columnls = load_csv(csv)

    for c in columnls:
        filtered_data, fst_list, pctcut = high_FST4comparison(data, c, pct)
        plotting_dist(fst_list, pctcut, csv + "_" + c)
        output_files(filtered_data, output_name + "_" + c.replace("Win_FST_", ""))

        print(c.replace("Win_FST_", ""), pctcut)








main_func("windowed_mean_5kbp_FST_averages_ChromX_genes.csv", 95)
main_func("windowed_mean_5kbp_FST_averages_Autosomes_genes.csv", 95)
main_func("windowed_max_5kbp_FST_averages_ChromX_genes.csv", 95)
main_func("windowed_max_5kbp_FST_averages_Autosomes_genes.csv", 95)







