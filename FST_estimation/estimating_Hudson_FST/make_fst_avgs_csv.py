import pandas as pd


def avgs_df(fstdf, avgcol):
    """
    Function takes a df with persite FST data for multiple comparisons.
    Negative FST values are clipped to zero. FST is averaged across comparisons.
    Returns a df with the columns CHROM POS avg_fst
    """

    #comparison columns
    cols = list(fstdf.columns)
    cols.remove("CHROM")
    cols.remove("POS")

    #clipping to zero
    fstdf.loc[:, cols] = fstdf.loc[:, cols].clip(lower=0)

    #averaging across comparisons
    fstdf[avgcol] = fstdf[cols].mean(axis=1)

    #dropping columns
    outdf = fstdf.drop(columns=cols)

    return outdf







#files to process
csvs = ["Zim_fst.csv", "ZS_fst.csv", "ZH_fst.csv", "ZC_fst.csv", "ZR_fst.csv", "Cos_fst.csv"]

#averaging autosomes
A_avgs = []
for csv in csvs:
    df = pd.read_csv(csv)
    df = df[df['CHROM'] != "X"]

    avgdf = avgs_df(df, "{pop}_FST_AVG".format(pop = csv.split("_")[0]))

    A_avgs.append(avgdf)


#averaging ChromX
X_avgs = []
for csv in csvs:
    df = pd.read_csv(csv)
    df = df[df['CHROM'] == "X"]

    avgdf = avgs_df(df, "{pop}_FST_AVG".format(pop = csv.split("_")[0]))

    X_avgs.append(avgdf)



#merging Autosome dfs
merged_A_df = A_avgs[0]
for df in A_avgs[1:]:
    merged_A_df = pd.merge(merged_A_df, df, on=["CHROM", "POS"], how="inner")


#merging ChromX dfs
merged_X_df = X_avgs[0]
for df in X_avgs[1:]:
    merged_X_df = pd.merge(merged_X_df, df, on=["CHROM", "POS"], how="inner")



merged_A_df.to_csv("FST_avergaes_Autosomes.csv", index=False)
merged_X_df.to_csv("FST_avergaes_ChromX.csv", index=False)