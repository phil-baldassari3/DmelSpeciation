import pandas as pd


def windowed_autosome_splitter(df):
    """
    splits an autosome df into searate chromosome dataframes
    """

    df_2L = df[df["Chrom"] == "2L"]
    df_2R = df[df["Chrom"] == "2R"]
    df_3L = df[df["Chrom"] == "3L"]
    df_3R = df[df["Chrom"] == "3R"]

    return df_2L, df_2R, df_3L, df_3R



def main_func(csvfile):

    windowed_autosomes = pd.read_csv(csvfile)

    data_2L, data_2R, data_3L, data_3R = windowed_autosome_splitter(windowed_autosomes)

    name2L = csvfile.replace("Autosomes", "Chrom2L")
    name2R = csvfile.replace("Autosomes", "Chrom2R")
    name3L = csvfile.replace("Autosomes", "Chrom3L")
    name3R = csvfile.replace("Autosomes", "Chrom3R")

    data_2L.to_csv(name2L, index=False)
    data_2R.to_csv(name2R, index=False)
    data_3L.to_csv(name3L, index=False)
    data_3R.to_csv(name3R, index=False)



main_func("windowed_mean_5kbp_FST_averages_Autosomes.csv")
main_func("windowed_max_5kbp_FST_averages_Autosomes.csv")