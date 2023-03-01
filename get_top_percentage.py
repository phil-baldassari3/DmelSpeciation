#importing modules
import os,sys
import pandas as pd
import numpy as np
from multiprocessing import Pool


###########Set Parameters################################################################################
#########################################################################################################

directory = "/Users/philipbaldassari/Desktop/annotated_windowed_fst"
output_directory = "/Users/philipbaldassari/Desktop/annotated_windowed_fst/top_percentage_csvs/"
top_percentage = 1

#########################################################################################################


def get_top_percent(infile, percentage, parameter, out_dir):
    """
    Function that takes a csv file with Fst data and extracts to top percentile inot a new dataframe.
    infile: {str} csv filename
    percentage: {int/float} percentage threshold you want to extact
    parameter: {str} name of column to sort by
    out_dir: {str} path to output directory
    """

    #reading csv to dataframe
    df = pd.read_csv(infile)

    #size of top percentage df
    top_len = int(len(df) * (percentage/100))

    #sorting df by parameter in descending order
    sorted_df = df.sort_values(by=parameter, ascending=False)

    #subset df
    top_df = sorted_df.iloc[0:top_len]

    #output file name
    name = "top{percent}_".format(percent=percentage) + infile

    #outputting file
    top_df.to_csv(out_dir + name, index=False)

    print("completed extraction of " + infile)



#########################################################################################################

os.chdir(directory)


#running the function
get_top_percent("windowed_1kbp_annotated_ZS_RAL_ZI_Fst_ChrX.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_ZS_RAL_ZI_Fst_autosomes.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_ZS_RAL_ZI_FR_SAfr_Fst_ChrX.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_ZS_RAL_ZI_FR_SAfr_Fst_autosomes.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_ZS_ZH_ZW_Fst_ChrX.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_ZS_ZH_ZW_Fst_autosomes.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_ZH_RAL_ZI_Fst_ChrX.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_ZH_RAL_ZI_Fst_autosomes.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_ZW_RAL_ZI_Fst_ChrX.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_ZW_RAL_ZI_Fst_autosomes.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_RAL_vs_FR_ZI_SAfr_Fst_ChrX.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_RAL_vs_FR_ZI_SAfr_Fst_autosomes.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_FR_vs_RAL_ZI_SAfr_Fst_ChrX.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_FR_vs_RAL_ZI_SAfr_Fst_autosomes.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_ZI_vs_RAL_FR_SAfr_Fst_ChrX.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_ZI_vs_RAL_FR_SAfr_Fst_autosomes.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_SAfr_vs_RAL_FR_ZI_Fst_ChrX.csv", top_percentage, "Avg_fst", output_directory)
get_top_percent("windowed_1kbp_annotated_SAfr_vs_RAL_FR_ZI_Fst_autosomes.csv", top_percentage, "Avg_fst", output_directory)

