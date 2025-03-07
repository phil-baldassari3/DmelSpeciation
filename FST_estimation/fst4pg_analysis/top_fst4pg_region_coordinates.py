import pandas as pd
import os


focalpop = "Zim"
comparison_pops = ["EAfr", "OOA.NW", "OOA.OW", "WCAfr", "SAfr", "Zam"]


for pop in comparison_pops:

    #loading csvs as dataframes
    df_X = pd.read_csv(f"{focalpop}/top_regions_{focalpop}_{pop}_ChrX.csv")

    df_2L = pd.read_csv(f"{focalpop}/top_regions_{focalpop}_{pop}_Chr2L.csv")
    df_2R = pd.read_csv(f"{focalpop}/top_regions_{focalpop}_{pop}_Chr2R.csv")
    df_3L = pd.read_csv(f"{focalpop}/top_regions_{focalpop}_{pop}_Chr3L.csv")
    df_3R = pd.read_csv(f"{focalpop}/top_regions_{focalpop}_{pop}_Chr3R.csv")


    #editing dfs
    df_X = df_X[["Start", "End"]]
    df_X.insert(0, "Chrom", "X")


    df_2L = df_2L[["Start", "End"]]
    df_2L.insert(0, "Chrom", "2L")
    
    df_2R = df_2R[["Start", "End"]]
    df_2R.insert(0, "Chrom", "2R")
    
    df_3L = df_3L[["Start", "End"]]
    df_3L.insert(0, "Chrom", "3L")
    
    df_3R = df_3R[["Start", "End"]]
    df_3R.insert(0, "Chrom", "3R")


    #concatenate autosomes
    df_autosomes = pd.concat([df_2L, df_2R, df_3L, df_3R], ignore_index=True)
    

    #saving dataframes
    df_X.to_csv(f"top_coordinates_{focalpop}_{pop}_ChrX.txt", index=False, header=False)
    df_autosomes.to_csv(f"top_coordinates_{focalpop}_{pop}_Autosomes.txt", index=False, header=False)


    #editing text files
    with open(f"top_coordinates_{focalpop}_{pop}_ChrX.txt", "r") as Xfile:
        Xstr = Xfile.read()
        Xstr = Xstr.replace("X,", "X:")
        Xstr = Xstr.replace(",", "..")
    with open(f"top_coordinates_{focalpop}_{pop}_ChrX.txt", "w") as Xfile:
        Xfile.write(Xstr)


    with open(f"top_coordinates_{focalpop}_{pop}_Autosomes.txt", "r") as Afile:
        Astr = Afile.read()
        Astr = Astr.replace("L,", "L:")
        Astr = Astr.replace("R,", "R:")
        Astr = Astr.replace(",", "..")
    with open(f"top_coordinates_{focalpop}_{pop}_Autosomes.txt", "w") as Afile:
        Afile.write(Astr)
    