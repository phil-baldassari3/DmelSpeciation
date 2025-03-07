import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def find_percentile_cutoffs(fst_ls, pct1, pct2):
    """
    Finds 2 percentile cutoffs from a list of FST values.
    pct1 and pct2 are teh chosen percentile cutoffs.
    75th percentile can be input as 75
    Returns 2 FST cutoff values
    """

    FST_cutoff1, FST_cutoff2 = np.percentile(fst_ls, [pct1, pct2])

    return FST_cutoff1, FST_cutoff2



def plot_cutoffs(fst_ls, cutoff1, cutoff2, title):
    """
    plots the FST distribution and the two percentile cutoff in standard adn log scale
    Saves 2 plots with filenames based on the given title
    """

    #number of bins, Freedman–Diaconis
    q25, q75 = np.percentile(fst_ls, [25, 75])
    bin_width = 2 * (q75 - q25) * len(fst_ls) ** (-1/3)
    bins = round((max(fst_ls) - min(fst_ls)) / bin_width)

    #plotting
    plt.figure(figsize=(10, 8))
    plt.hist(fst_ls, density=True, bins=bins)
    plt.axvline(x=cutoff1, color='r', lw=1, linestyle='--')
    plt.axvline(x=cutoff2, color='r', lw=1, linestyle='--')
    plt.axvspan(cutoff1, cutoff2, alpha=0.3, color='gray')
    plt.ylabel('Density')
    plt.xlabel("FST")
    plt.title(title)

    plt.savefig(f"{title}_dist.png")

    plt.clf()

    #plotting log
    plt.figure(figsize=(10, 8))
    plt.hist(fst_ls, density=True, bins=bins)
    plt.axvline(x=cutoff1, color='r', lw=1, linestyle='--')
    plt.axvline(x=cutoff2, color='r', lw=1, linestyle='--')
    plt.axvspan(cutoff1, cutoff2, alpha=0.3, color='gray')
    plt.yscale('log')
    plt.ylabel('Density')
    plt.xlabel("FST")
    plt.title(title)

    plt.savefig(f"{title}_dist_logscale.png")

    plt.clf()



def make_diff_Z_df(df, focalpop, cospop, cutoff_cos, cutoff_focal):
    """
    Finds significantly differentiated Zimbabwe SNPs by
    filtering the dataframe for Zim FST above the focal cutoff and Cos FST below the cosmopolitan cutoff
    Returns the filtered dataframe
    """

    #remove unused columns
    cols = list(df.columns)
    cols.remove("CHROM")
    cols.remove("POS")
    cols.remove(f"{focalpop}_FST_AVG")
    cols.remove(f"{cospop}_FST_AVG")

    df = df.drop(columns=cols)

    #filtering df
    filtered_df = df[(df[f"{focalpop}_FST_AVG"] > cutoff_focal) & (df[f"{cospop}_FST_AVG"] < cutoff_cos)]

    return filtered_df



def main_func(csvfile, focal_pop, cos_pop, plottitle):
    """
    Takes in a FST averages csv, a focal and cosmopolitan population and outputs a csv of SNPs
    that are above the 95th percentile of FST values for the focal population and below the 
    75th percentile of FST for the cosmopolitan population. The fucntion also plots the FST
    distribution with the 75th and 95th percentile cutoffs
    """
    #opening data
    fst_data = pd.read_csv(csvfile)

    #making total list of FST values
    focal_ls = fst_data[f"{focal_pop}_FST_AVG"].to_list()
    cos_ls = fst_data[f"{cos_pop}_FST_AVG"].to_list()
    total_ls = focal_ls + cos_ls

    #finding percentile cutoffs
    cos_cutoff, focal_cutoff = find_percentile_cutoffs(total_ls, 75, 95)
    print(f"{cos_pop}_FST < {cos_cutoff} and {focal_pop}_FST > {focal_cutoff}")

    #plotting distibution
    plot_cutoffs(total_ls, cos_cutoff, focal_cutoff, plottitle)

    #filtering data with cutoffs
    filtered_SNPs = make_diff_Z_df(fst_data, focal_pop, cos_pop, cos_cutoff, focal_cutoff)

    #saving df
    if "ChromX" in csvfile:
        outputname = f"{focal_pop}_diff_SNPs_ChromX.csv"
    elif "Autosomes" in csvfile:
        outputname = f"{focal_pop}_diff_SNPs_Autosomes.csv"
    else:
        print("cannot parse name")

    filtered_SNPs.to_csv(outputname, index=False)

    










print("ChromX")
main_func("FST_avergaes_ChromX.csv", "Zim", "Cos", "Zim and Cos ChromX FST Averages with 75th and 95th Percentiles")
main_func("FST_avergaes_ChromX.csv", "ZS", "Cos", "ZS and Cos ChromX FST Averages with 75th and 95th Percentiles")
main_func("FST_avergaes_ChromX.csv", "ZH", "Cos", "ZH and Cos ChromX FST Averages with 75th and 95th Percentiles")
main_func("FST_avergaes_ChromX.csv", "ZC", "Cos", "ZC and Cos ChromX FST Averages with 75th and 95th Percentiles")
main_func("FST_avergaes_ChromX.csv", "ZR", "Cos", "ZR and Cos ChromX FST Averages with 75th and 95th Percentiles")
print("Autosomes")
main_func("FST_avergaes_Autosomes.csv", "Zim", "Cos", "Zim and Cos Autosomal FST Averages with 75th and 95th Percentiles")
main_func("FST_avergaes_Autosomes.csv", "ZS", "Cos", "ZS and Cos Autosomal FST Averages with 75th and 95th Percentiles")
main_func("FST_avergaes_Autosomes.csv", "ZH", "Cos", "ZH and Cos Autosomal FST Averages with 75th and 95th Percentiles")
main_func("FST_avergaes_Autosomes.csv", "ZC", "Cos", "ZC and Cos Autosomal FST Averages with 75th and 95th Percentiles")
main_func("FST_avergaes_Autosomes.csv", "ZR", "Cos", "ZR and Cos Autosomal FST Averages with 75th and 95th Percentiles")





