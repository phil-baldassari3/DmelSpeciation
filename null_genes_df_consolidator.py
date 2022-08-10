#script consolidates the null condition Fst dfs with the gene hits dfs and saves list to genes to run through GO

import os,sys
import pandas as pd
import numpy as np
from multiprocessing import Pool

#setting directory
directory = "/Users/philipbaldassari/Desktop/zim-cos_downsampled/null_condition"
os.chdir(directory)

'''
#get list of txt list files for parallel pool mapping
file_list = []

for file in os.listdir(directory):
    if file.endswith('.csv') and 'null' in str(file) and 'gene' not in str(file):
        file_list.append(file)
    else:
        continue

print(file_list)
'''

#df consolidator function
def df_consolidator(infile):

    #opening sites df
    site_df = pd.read_csv(infile)

    for file in os.listdir(directory):
        if file.endswith('.csv') and str(infile) in str(file) and 'gene_hits' in str(file):

            #opening coresponding gene df
            gene_df = pd.read_csv(file)

            print("corresponding dataframes opened for ", infile, " now merging...")

            #merging dfs
            consolidated_df = pd.merge(site_df, gene_df, left_index=True, right_index=True)

            #getting list of genes and cleaning
            gene_list = gene_df['FBgn'].tolist()
            gene_list4GO = [value for value in gene_list if value != "None"]

            #saving gene list to file
            with open('gene4GO_{name}.txt'.format(name=infile.split('.')[0]), 'w') as gene_ls:
                gene_ls.write(', '.join(gene_list4GO))

            print("saved gene list file for GO for ", infile)

            #saving consolidated df to file
            consolidated_df.to_csv('sites_genes_{name}.csv'.format(name=infile.split('.')[0]), index=False)
            
            print("saved consolidated dataframe for ", infile)

            break

        else:
            continue





#for parallel mapping
csv_list = ['null_Fst_ZS_ZH_ZW_ChrX.csv', 'null_Fst_Zim_RAL_ZI_Chr2R.csv', 'null_Fst_Zim_RAL_ZI_autosome.csv', 'null_Fst_ZS_RAL_ZI_FR_SAfr_autosome.csv', 'null_Fst_ZH_RAL_ZI_ChrX.csv', 'null_Fst_ZS_ZH_ZW_Chr3R.csv', 'null_Fst_ZS_ZH_ZW_Chr2R.csv', 'null_Fst_Zim_RAL_ZI_Chr3R.csv', 'null_Fst_ZS_RAL_ZI_Chr3R.csv', 'null_Fst_ZW_RAL_ZI_Chr2R.csv', 'null_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L.csv', 'null_Fst_ZH_RAL_ZI_Chr2R.csv', 'null_Fst_ZS_RAL_ZI_ChrX.csv', 'null_Fst_ZH_RAL_ZI_Chr3R.csv', 'null_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L.csv', 'null_Fst_ZW_RAL_ZI_Chr3R.csv', 'null_Fst_ZS_RAL_ZI_FR_SAfr_ChrX.csv', 'null_Fst_ZS_RAL_ZI_Chr2R.csv', 'null_Fst_ZS_RAL_ZI_Chr3L.csv', 'null_Fst_Zim_RAL_ZI_ChrX.csv', 'null_Fst_ZW_RAL_ZI_Chr2L.csv', 'null_Fst_ZW_RAL_ZI_autosome.csv', 'null_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R.csv', 'null_Fst_ZH_RAL_ZI_Chr2L.csv', 'null_Fst_ZH_RAL_ZI_Chr3L.csv', 'null_Fst_ZS_RAL_ZI_autosome.csv', 'null_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R.csv', 'null_Fst_ZW_RAL_ZI_Chr3L.csv', 'null_Fst_ZS_RAL_ZI_Chr2L.csv', 'null_Fst_Zim_RAL_ZI_Chr2L.csv', 'null_Fst_ZS_ZH_ZW_autosome.csv', 'null_Fst_ZS_ZH_ZW_Chr3L.csv', 'null_Fst_ZS_ZH_ZW_Chr2L.csv', 'null_Fst_ZH_RAL_ZI_autosome.csv', 'null_Fst_Zim_RAL_ZI_Chr3L.csv', 'null_Fst_ZW_RAL_ZI_ChrX.csv']

#run in parallel
def run_in_parallel():
    pool = Pool(processes=8)
    pool.map(df_consolidator, csv_list)


if __name__ == '__main__':
    run_in_parallel()


       

            
