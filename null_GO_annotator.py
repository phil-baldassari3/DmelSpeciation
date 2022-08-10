#script to annotate the null genes for neurogenesis and mating behavior GO terms and child terms
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
    if file.endswith('.csv') and 'sites_genes' in str(file):
        file_list.append(file)
    else:
        continue

print(file_list)
'''

###Flybase GO lists###
#GO neurogenesis file
FB_neurogenesis_file = open('flybase_GO_neurogenesis_genes.txt', "r")
FB_neurogenesis_genes = FB_neurogenesis_file.read()
#to list
FB_neurogenesis_GO = FB_neurogenesis_genes.split("\n")
FB_neurogenesis_file.close()

#GO mating behavior file
FB_mating_behavior_file = open('flybase_GO_matingbehavior_genes.txt', "r")
FB_mating_behavior_genes = FB_mating_behavior_file.read()
#to list
FB_mating_behavior_GO = FB_mating_behavior_genes.split("\n")
FB_mating_behavior_file.close()

#GO male mating behavior file
FB_male_mating_behavior_file = open('flybase_GO_malemating_genes.txt', "r")
FB_male_mating_behavior_genes = FB_male_mating_behavior_file.read()
#to list
FB_male_mating_behavior_GO = FB_male_mating_behavior_genes.split("\n")
FB_male_mating_behavior_file.close()

#GO female mating behavior file
FB_female_mating_behavior_file = open('flybase_GO_femalemating_genes.txt', "r")
FB_female_mating_behavior_genes = FB_female_mating_behavior_file.read()
#to list
FB_female_mating_behavior_GO = FB_female_mating_behavior_genes.split("\n")
FB_female_mating_behavior_file.close()



###Flymine GO lists###
#GO neurogenesis file
FM_neurogenesis_df = pd.read_csv('flymine_GO_neurogenesis_genes.csv')
#to list
FM_neurogenesis_GO = FM_neurogenesis_df['FBgn'].tolist()

#GO mating behavior file
FM_mating_behavior_df = pd.read_csv('flymine_GO_matingbehavior_genes.csv')
#to list
FM_mating_behavior_GO = FM_mating_behavior_df['FBgn'].tolist()

#GO male mating behavior file
FM_male_mating_behavior_df = pd.read_csv('flymine_GO_malemating_genes.csv')
#to list
FM_male_mating_behavior_GO = FM_male_mating_behavior_df['FBgn'].tolist()

#GO female mating behavior file
FM_female_mating_behavior_df = pd.read_csv('flymine_GO_femalemating_genes.csv')
#to list
FM_female_mating_behavior_GO = FM_female_mating_behavior_df['FBgn'].tolist()

print("Read GO files to lists")

wait = input("Press Enter to continue.")


#GO annotator
def GO_annotator(infile):

    #openning df
    df = pd.read_csv(infile)

    #gene list
    gene_list = df['FBgn'].tolist()

    print("opened df from ", infile, " and read genes to list")

    #setting lists
    GO_neurogenesis_FB = []
    GO_mating_behavior_FB = []
    GO_male_mating_behavior_FB = []
    GO_female_mating_behavior_FB = []
    GO_neurogenesis_FM = []
    GO_mating_behavior_FM = []
    GO_male_mating_behavior_FM = []
    GO_female_mating_behavior_FM = []

    print("looping through the gene list from ", infile)

    #looping through null genes
    for gene in gene_list:
        
        if gene in FB_neurogenesis_GO:
            GO_neurogenesis_FB.append('Neurogenesis')
        else:
            GO_neurogenesis_FB.append('')
        

        if gene in FB_mating_behavior_GO:
            GO_mating_behavior_FB.append('Mating_behavior')
        else:
            GO_mating_behavior_FB.append('')
        

        if gene in FB_male_mating_behavior_GO:
            GO_male_mating_behavior_FB.append('Male_mating_behavior')
        else:
            GO_male_mating_behavior_FB.append('')
        

        if gene in FB_female_mating_behavior_GO:
            GO_female_mating_behavior_FB.append('Female_mating_behavior')
        else:
            GO_female_mating_behavior_FB.append('')
        



        if gene in FM_neurogenesis_GO:
            GO_neurogenesis_FM.append('Neurogenesis')
        else:
            GO_neurogenesis_FM.append('')
        

        if gene in FM_mating_behavior_GO:
            GO_mating_behavior_FM.append('Mating_behavior')
        else:
            GO_mating_behavior_FM.append('')
        

        if gene in FM_male_mating_behavior_GO:
            GO_male_mating_behavior_FM.append('Male_mating_behavior')
        else:
            GO_male_mating_behavior_FM.append('')
        

        if gene in FM_female_mating_behavior_GO:
            GO_female_mating_behavior_FM.append('Female_mating_behavior')
        else:
            GO_female_mating_behavior_FM.append('')

    #making dataframe
    dictionary = {'FB_GO:0022008' : GO_neurogenesis_FB, 'FM_GO:0022008' : GO_neurogenesis_FM, 'FB_GO:0007617' : GO_mating_behavior_FB, 'FM_GO:0007617' : GO_mating_behavior_FM, 'FB_GO:0060179' : GO_male_mating_behavior_FB, 'FM_GO:0060179' : GO_male_mating_behavior_FM, 'FB_GO:0060180' : GO_female_mating_behavior_FB, 'FM_GO:0060180' : GO_female_mating_behavior_FM}

    #making df
    GO_df = pd.DataFrame(dictionary)

    print("GO df made from ", infile)

    #merging dataframes
    merged_df = pd.merge(df, GO_df, left_index=True, right_index=True)

    print("merging dataframes")

    #saving dataframe as csv
    merged_df.to_csv('GO_{file}.csv'.format(file=infile), index=False)

    print("saved GO_{file}.csv".format(file=infile))





#for parallel mapping
csv_list = ['sites_genes_null_Fst_ZS_ZH_ZW_ChrX.csv', 'sites_genes_null_Fst_ZS_RAL_ZI_ChrX.csv', 'sites_genes_null_Fst_Zim_RAL_ZI_Chr2L.csv', 'sites_genes_null_Fst_Zim_RAL_ZI_Chr3L.csv', 'sites_genes_null_Fst_ZW_RAL_ZI_Chr3L.csv', 'sites_genes_null_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R.csv', 'sites_genes_null_Fst_ZS_ZH_ZW_autosome.csv', 'sites_genes_null_Fst_ZS_RAL_ZI_Chr2L.csv', 'sites_genes_null_Fst_Zim_RAL_ZI_autosome.csv', 'sites_genes_null_Fst_ZS_ZH_ZW_Chr2R.csv', 'sites_genes_null_Fst_ZH_RAL_ZI_Chr3L.csv', 'sites_genes_null_Fst_ZH_RAL_ZI_ChrX.csv', 'sites_genes_null_Fst_ZH_RAL_ZI_Chr2L.csv', 'sites_genes_null_Fst_ZS_ZH_ZW_Chr3R.csv', 'sites_genes_null_Fst_ZS_RAL_ZI_Chr3L.csv', 'sites_genes_null_Fst_Zim_RAL_ZI_ChrX.csv', 'sites_genes_null_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R.csv', 'sites_genes_null_Fst_ZS_RAL_ZI_FR_SAfr_autosome.csv', 'sites_genes_null_Fst_ZW_RAL_ZI_Chr2L.csv', 'sites_genes_null_Fst_ZW_RAL_ZI_Chr3R.csv', 'sites_genes_null_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L.csv', 'sites_genes_null_Fst_ZS_RAL_ZI_Chr2R.csv', 'sites_genes_null_Fst_ZH_RAL_ZI_autosome.csv', 'sites_genes_null_Fst_ZW_RAL_ZI_ChrX.csv', 'sites_genes_null_Fst_ZS_ZH_ZW_Chr2L.csv', 'sites_genes_null_Fst_ZH_RAL_ZI_Chr3R.csv', 'sites_genes_null_Fst_ZS_RAL_ZI_FR_SAfr_ChrX.csv', 'sites_genes_null_Fst_ZH_RAL_ZI_Chr2R.csv', 'sites_genes_null_Fst_ZS_ZH_ZW_Chr3L.csv', 'sites_genes_null_Fst_ZS_RAL_ZI_Chr3R.csv', 'sites_genes_null_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L.csv', 'sites_genes_null_Fst_ZW_RAL_ZI_Chr2R.csv', 'sites_genes_null_Fst_ZS_RAL_ZI_autosome.csv', 'sites_genes_null_Fst_Zim_RAL_ZI_Chr2R.csv', 'sites_genes_null_Fst_Zim_RAL_ZI_Chr3R.csv', 'sites_genes_null_Fst_ZW_RAL_ZI_autosome.csv']
#run in parallel
def run_in_parallel():
    pool = Pool(processes=8)
    pool.map(GO_annotator, csv_list)


if __name__ == '__main__':
    run_in_parallel()


        

