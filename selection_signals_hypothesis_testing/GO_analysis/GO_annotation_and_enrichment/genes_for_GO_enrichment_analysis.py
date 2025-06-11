import pandas as pd
import numpy as np


def percentage_genes_from_windows(csvfile, percentage, direction):
    """
    Function takes in a csv file with windowed genome statistics data with columns:
    CHROM, START, Statistics_data..., FBgn, gene_symbol
    and saves lists of genes represented by the top or bottom percentage of windows,
    one list of each statistics column.

    The percentage is computed as either the top P% or bottom P% for each chromosome
    separately to account of linkage, heritability, and chromosomal differences.

    Files are saved as text files, filenames generated from input csv name and column
    names.

    Arguments:
    csvfile (str): input csv file
    percentage (int, float): percentage of windows to extract genes from
    direction (str): "top" or "bottom". the direction of the percentage. extracts either the "top" P% or "bottom" P% of windows
    """

    print("\n")

    #loading data
    df = pd.read_csv(csvfile, keep_default_na=False)

    #finding percentile threshold
    if direction == "top":
        percent = 100 - percentage
    elif direction == "bottom":
        percent = percentage
    else:
        raise ValueError("Wrong direction set. Use either 'top' or 'bottom'.")

    #getting chromsome names
    chroms = list(df["CHROM"].unique())

    #finding data column names
    annotation_cols = ["CHROM", "START", "END", "FBgn", "gene_symbol"]
    data_cols = list(df.columns)
    data_cols = [i for i in data_cols if i not in annotation_cols]

    #looping through data columns
    for colname in data_cols:

        genes = []
        
        print(colname)

        #looping through chrom dfs
        for chrom in chroms:

            #extracting chrom
            filtered_df = df[df["CHROM"] == chrom]

            #finding percentile
            stat_values = filtered_df[colname].to_list()
            pct = np.percentile(stat_values, percent)

            #filtering stat values
            if direction == "top":
                filtered_df = filtered_df[filtered_df[colname] >= pct]
            elif direction == "bottom":
                filtered_df = filtered_df[filtered_df[colname] <= pct]

            #extract genes to list
            tempgenes = filtered_df["FBgn"].to_list()
            tempgenes = [i for i in tempgenes if i != "None"]
            tempgenes = "; ".join(tempgenes)
            chrom_genes = list(set(tempgenes.split("; ")))
            
            #appending genes to master list
            genes += chrom_genes

            print(f"Chrom {chrom} Threshold: {pct}")

        print(f"Number of genes: {len(genes)}")

        #saving genes file for this column
        outputname = csvfile.split("/")[-1]
        outputname = outputname.replace(".csv", ".txt")
        outputname = outputname.replace("selected_pops_", "")
        outputname = f"{colname}_{direction}_{percentage}%_genes_" + outputname

        with open(outputname, "w") as outfile:
            outfile.write("\n".join(genes))




def percentage_genes_from_per_gene(csvfile, percentage, direction):
    """
    Function takes in a csv file with per-gene genome statistics data with columns:
    FBgn, gene_symbol, Statistics_data...
    and saves lists of genes represented by the top or bottom percentage,
    one list of each statistics column. The percentage is computed as either the top
    P% or bottom P%.

    Files are saved as text files, filenames generated from input csv name and column
    names.

    Arguments:
    csvfile (str): input csv file
    percentage (int, float): percentage of genes to extract
    direction (str): "top" or "bottom". the direction of the percentage. extracts either the "top" P% or "bottom" P% of genes
    """

    print("\n")

    #loading data
    df = pd.read_csv(csvfile, keep_default_na=False)

    #finding percentile threshold
    if direction == "top":
        percent = 100 - percentage
    elif direction == "bottom":
        percent = percentage
    else:
        raise ValueError("Wrong direction set. Use either 'top' or 'bottom'.")

    #finding data column names
    annotation_cols = ["FBgn", "gene_symbol"]
    data_cols = list(df.columns)
    data_cols = [i for i in data_cols if i not in annotation_cols]

    #looping through data columns
    for colname in data_cols:

        genes = []
        
        print(colname)

        #finding percentile
        stat_values = df[colname].to_list()
        pct = np.percentile(stat_values, percent)

        #filtering stat values
        if direction == "top":
            filtered_df = df[df[colname] >= pct]
        elif direction == "bottom":
            filtered_df = df[df[colname] <= pct]

        #extract genes to list
        genes = filtered_df["FBgn"].to_list()

        print(f"Threshold: {pct}")

        print(f"Number of genes: {len(genes)}")

        #saving genes file for this column
        outputname = csvfile.split("/")[-1]
        outputname = outputname.replace(".csv", ".txt")
        outputname = outputname.replace("selected_pops_", "")
        outputname = f"{colname}_{direction}_{percentage}%_genes_" + outputname

        with open(outputname, "w") as outfile:
            outfile.write("\n".join(genes))





#per-gene FST
percentage_genes_from_per_gene("FST/FST_data/selected_pops_per_gene_fst.csv", 10, "top")
percentage_genes_from_per_gene("FST/FST_data/selected_pops_per_gene_fst_NormByGeneLength.csv", 10, "top")

#windowed FST
percentage_genes_from_windows("FST/FST_data/windowed_mean_10kbp_selected_pops_fst_mapped2genes.csv", 10, "top")


#per-gene PI
percentage_genes_from_per_gene("Pi/Pi_data/selected_pops_per_gene_pi.csv", 10, "bottom")


#windwoed PI
percentage_genes_from_windows("Pi/Pi_data/selected_pops_windowed_10kbp_pi_mapped2genes.csv", 10, "bottom")


#theta
percentage_genes_from_windows("Theta/Theta_data/selected_pops_windowed_10kbp_theta_mapped2genes.csv", 10, "bottom")


#tajima's D
percentage_genes_from_windows("TajimasD/TajimasD_data/selected_pops_windowed_10kbp_TajimasD_mapped2genes.csv", 10, "bottom")






