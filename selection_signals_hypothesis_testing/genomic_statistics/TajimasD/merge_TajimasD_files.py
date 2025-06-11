import pandas as pd


def load_tajima_data(pop):

    df = pd.read_csv(f"{pop}_TajimasD.tsv", sep="\t")

    df.drop(["N_Sites", "N_SNPs"], axis=1, inplace=True)

    df.rename(columns={"BIN_START":"START", "BIN_END":"END", "TajimaD":f"{pop}_TajimaD"}, inplace=True)

    df["START"] = df["START"] + 1

    return df


ZSdf = load_tajima_data("ZS")
ZHdf = load_tajima_data("ZH")
RALdf = load_tajima_data("RAL")
FRdf = load_tajima_data("FR")
ZIdf = load_tajima_data("ZI")

dfs = [ZSdf, ZHdf, RALdf, FRdf, ZIdf]
merged_df = dfs[0]
for df in dfs[1:]:
    merged_df = pd.merge(merged_df, df, on=["CHROM", "START", "END"], how="inner")


merged_df.to_csv("selected_pops_windowed_10kbp_TajimasD.csv", index=False)