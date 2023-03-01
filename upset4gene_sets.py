import os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import upsetplot



def plot_upset_from_gene_csvs(csv_ls, plotname, subset="all"):
    """
    makes upset plot from lists of genes from the top gene csv files. A list of csv files is the first argument. The files will be opened and the gene lists extracted.
    Pass the title of the plot you wish to output as plotname
    subset is an optional argument if you would like to subset the gene lists by a specific paramenter of the dataframe. Pass the exact parameter name to subset
    Also returns a modified upset dataframe
    """

    gene_ls_dict = {}

    for i in csv_ls:

        #name of the set
        comparison_name = '_'.join(i.split('_')[4:]).replace('_Fst_ChrX.csv', '')
        comparison_name = comparison_name.replace('_Fst_autosomes.csv', '')

        #opening and subseting df if necessary
        if subset == "all":
            df = pd.read_csv(i)
        else:
            full_df = pd.read_csv(i)
            df = full_df[full_df.apply(lambda row: row.astype(str).str.match(subset, case=False).any(), axis=1)]

        #updating dictionary
        gene_ls_dict.update({comparison_name: df["FBgn"].to_list()})


    #converting dictionary to upset format
    gene_sets = upsetplot.from_contents(gene_ls_dict)

    #plotting
    p = upsetplot.UpSet(gene_sets, subset_size='count', sort_by='-degree', sort_categories_by='-input').plot()
#sort_by='cardinality'
    plt.title(plotname)
    plt.savefig(plotname + '_upset.png')


    return gene_sets.reset_index()
    

        
def extract_subset(df, bool_ls):
    """
    Function takes a upset dataframe and a bool list defining the rows of the dataframe to select and returns a list of items of the seelcted subset
    Note that the bool list needs to be in the same order as the catagories passed to plot_upset_from_gene_csvs()
    """
    #find how many rows to subset from
    length = len(bool_ls)

    #getting subset
    subset_df = df.loc[(df.iloc[:, :length] == bool_ls).all(axis=1)]

    #getting list of genes
    subset_genes = subset_df['id'].to_list()

    return subset_genes


def extract_all_subsets(df, num_of_catagories):
    """Function takes in the upset dataframe and the number of categorites and outputs a dictionary of lists of all subsets"""

    #unique list of booleans that describe the identiy of each item in each catagory
    unique_bools = df.iloc[:, :num_of_catagories].drop_duplicates().values.tolist()

    subsets_dict = {}
    for i in unique_bools:
        subset = extract_subset(df, i)

        subsets_dict.update({str(i): subset})

    return subsets_dict








chrx = ['top1_windowed_10kbp_genes_ZS_vs_RAL_ZI_Fst_ChrX.csv', 'top1_windowed_10kbp_genes_ZS_vs_RAL_ZI_FR_SAfr_Fst_ChrX.csv', 'top1_windowed_10kbp_genes_ZS_vs_ZH_ZW_Fst_ChrX.csv', 'top1_windowed_10kbp_genes_ZH_vs_RAL_ZI_Fst_ChrX.csv', 'top1_windowed_10kbp_genes_ZW_vs_RAL_ZI_Fst_ChrX.csv', 'top1_windowed_10kbp_genes_RAL_vs_FR_ZI_SAfr_Fst_ChrX.csv', 'top1_windowed_10kbp_genes_FR_vs_RAL_ZI_SAfr_Fst_ChrX.csv', 'top1_windowed_10kbp_genes_ZI_vs_RAL_FR_SAfr_Fst_ChrX.csv', 'top1_windowed_10kbp_genes_SAfr_vs_RAL_FR_ZI_Fst_ChrX.csv']
autosomes = ['top1_windowed_10kbp_genes_ZS_vs_RAL_ZI_Fst_autosomes.csv', 'top1_windowed_10kbp_genes_ZS_vs_RAL_ZI_FR_SAfr_Fst_autosomes.csv', 'top1_windowed_10kbp_genes_ZS_vs_ZH_ZW_Fst_autosomes.csv', 'top1_windowed_10kbp_genes_ZH_vs_RAL_ZI_Fst_autosomes.csv', 'top1_windowed_10kbp_genes_ZW_vs_RAL_ZI_Fst_autosomes.csv', 'top1_windowed_10kbp_genes_RAL_vs_FR_ZI_SAfr_Fst_autosomes.csv', 'top1_windowed_10kbp_genes_FR_vs_RAL_ZI_SAfr_Fst_autosomes.csv', 'top1_windowed_10kbp_genes_ZI_vs_RAL_FR_SAfr_Fst_autosomes.csv', 'top1_windowed_10kbp_genes_SAfr_vs_RAL_FR_ZI_Fst_autosomes.csv']


df_chrx_win10 = plot_upset_from_gene_csvs(chrx, "Top_Genes_ChrX_windowed10kbp")
df_chrx_neuro_win10 = plot_upset_from_gene_csvs(chrx, "Top_Neurogenesis_Genes_ChrX_windowed10kbp", subset='Neurogenesis')
df_autosome_win10 = plot_upset_from_gene_csvs(autosomes, "Top_Genes_Autosomes_windowed10kbp")
df_autosome_neuro_win10 = plot_upset_from_gene_csvs(autosomes, "Top_Neurogenesis_Genes_Autosomes_windowed10kbp", subset='Neurogenesis')






chrx = ['top1_windowed_1kbp_genes_ZS_vs_RAL_ZI_Fst_ChrX.csv', 'top1_windowed_1kbp_genes_ZS_vs_RAL_ZI_FR_SAfr_Fst_ChrX.csv', 'top1_windowed_1kbp_genes_ZS_vs_ZH_ZW_Fst_ChrX.csv', 'top1_windowed_1kbp_genes_ZH_vs_RAL_ZI_Fst_ChrX.csv', 'top1_windowed_1kbp_genes_ZW_vs_RAL_ZI_Fst_ChrX.csv', 'top1_windowed_1kbp_genes_RAL_vs_FR_ZI_SAfr_Fst_ChrX.csv', 'top1_windowed_1kbp_genes_FR_vs_RAL_ZI_SAfr_Fst_ChrX.csv', 'top1_windowed_1kbp_genes_ZI_vs_RAL_FR_SAfr_Fst_ChrX.csv', 'top1_windowed_1kbp_genes_SAfr_vs_RAL_FR_ZI_Fst_ChrX.csv']
autosomes = ['top1_windowed_1kbp_genes_ZS_vs_RAL_ZI_Fst_autosomes.csv', 'top1_windowed_1kbp_genes_ZS_vs_RAL_ZI_FR_SAfr_Fst_autosomes.csv', 'top1_windowed_1kbp_genes_ZS_vs_ZH_ZW_Fst_autosomes.csv', 'top1_windowed_1kbp_genes_ZH_vs_RAL_ZI_Fst_autosomes.csv', 'top1_windowed_1kbp_genes_ZW_vs_RAL_ZI_Fst_autosomes.csv', 'top1_windowed_1kbp_genes_RAL_vs_FR_ZI_SAfr_Fst_autosomes.csv', 'top1_windowed_1kbp_genes_FR_vs_RAL_ZI_SAfr_Fst_autosomes.csv', 'top1_windowed_1kbp_genes_ZI_vs_RAL_FR_SAfr_Fst_autosomes.csv', 'top1_windowed_1kbp_genes_SAfr_vs_RAL_FR_ZI_Fst_autosomes.csv']


df_chrx_win1 = plot_upset_from_gene_csvs(chrx, "Top_Genes_ChrX_windowed1kbp")
df_chrx_neuro_win1 = plot_upset_from_gene_csvs(chrx, "Top_Neurogenesis_Genes_ChrX_windowed1kbp", subset='Neurogenesis')
df_autosome_win1 = plot_upset_from_gene_csvs(autosomes, "Top_Genes_Autosomes_windowed1kbp")
df_autosome_neuro_win1 = plot_upset_from_gene_csvs(autosomes, "Top_Neurogenesis_Genes_Autosomes_windowed1kbp", subset='Neurogenesis')




bl = [True, True, True, False, False, False, False, False, False]
l = extract_subset(df_chrx_win10, bl)
print(l)
print(len(l))



d = extract_all_subsets(df_chrx_win10, 9)
df_test = (pd.DataFrame.from_dict(d, orient='index')).transpose()

df_test.to_csv('test.csv', index=False)
