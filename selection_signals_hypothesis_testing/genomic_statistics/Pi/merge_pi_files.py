import pandas as pd
import os

def load_pi_data(file, endswith):

    df = pd.read_csv(file, sep="\t")

    #extracting population name
    popname = file.split(".")[0]

    #rename pi column
    df.rename(columns={"PI": f"{popname}_PI"}, inplace=True)

    #extra steps for windowed data
    if endswith == ".windowed.pi":
        df.rename(columns={"BIN_START": "START", "BIN_END": "END"}, inplace=True)
        df.drop("N_VARIANTS", axis=1, inplace=True)

    return df


def merge_pi_data(endswith):

    #listing dataframes
    dfs = []

    for file in os.listdir():
        if file.endswith(endswith):
            
            pidf = load_pi_data(file, endswith)
            dfs.append(pidf)


    #starting df
    merged_df = dfs[0]

    if endswith == ".sites.pi":
        #looping and merging remaining dfs
        for df in dfs[1:]:
            merged_df = pd.merge(merged_df, df, on=["CHROM", "POS"], how="inner")

    elif endswith == ".windowed.pi":
        #looping and merging remaining dfs
        for df in dfs[1:]:
            merged_df = pd.merge(merged_df, df, on=["CHROM", "START", "END"], how="outer")
        merged_df.fillna(0, inplace=True)
        merged_df.sort_values(by=["CHROM", "START"], inplace=True)

    else:
        print("oops you did something wrong")

    return merged_df





persite_merged = merge_pi_data(".sites.pi")
windowed_merged = merge_pi_data(".windowed.pi")

persite_merged.to_csv("selected_pops_sites_pi.csv", index=False)
windowed_merged.to_csv("selected_pops_windowed_10kbp_pi.csv", index=False)