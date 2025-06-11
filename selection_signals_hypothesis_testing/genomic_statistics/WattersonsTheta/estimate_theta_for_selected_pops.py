import pandas as pd
import numpy as np


def estimate_theta(pop, num_samples):

    #load snpden data
    df = pd.read_csv(f"{pop}.snpden", sep="\t")

    #compute harmonic number a
    diploid_sample_size = num_samples * 2
    summation_range = list(range(1, diploid_sample_size))
    harmonics = [1/i for i in summation_range]
    a = sum(harmonics)

    #editing starts and adding stops
    df.rename(columns={"BIN_START": "START"}, inplace=True)
    df["END"] = df["START"] + 10000
    df["START"] = df["START"] + 1

    #compute theta
    df[f"{pop}_Theta"] = df["SNP_COUNT"] / a

    #drop unnecessary columns
    df.drop("VARIANTS/KB", axis=1, inplace=True)
    df.drop("SNP_COUNT", axis=1, inplace=True)

    return df



#estimate theta
ZS_df = estimate_theta("ZS", 13)
ZH_df = estimate_theta("ZH", 14)
RAL_df = estimate_theta("RAL", 84)
FR_df = estimate_theta("FR", 40)
ZI_df = estimate_theta("ZI", 51)


dfs = [ZS_df, ZH_df, RAL_df, FR_df, ZI_df]
merged_df = dfs[0]
for df in dfs[1:]:
    merged_df = pd.merge(merged_df, df, on=["CHROM", "START", "END"], how="outer")


merged_df.to_csv("selected_pops_windowed_10kbp_theta.csv", index=False)