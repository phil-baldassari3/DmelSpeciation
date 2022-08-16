#Script to permutate over the null set to count relevantly functional genes

#importing modules
import os,sys
import pandas as pd
import numpy as np
from multiprocessing import Pool
import random

#setting directory
directory = "/Users/philipbaldassari/Desktop/zim-cos_downsampled/null_condition"
os.chdir(directory)

'''
#get list of txt list files for parallel pool mapping
file_list = []

for file in os.listdir(directory):
    if file.endswith('.csv') and 'GO_sites_genes' in str(file):
        file_list.append(file)
    else:
        continue

print(file_list)
'''

#1% sample sizes
site_sample_size = {'ZS_RAL_ZI_Chr2L':15503, 'ZH_RAL_ZI_autosome':12211, 'ZW_RAL_ZI_Chr3L':7, 'ZS_RAL_ZI_FR_SAfr_Chr3L':14902, 'ZS_ZH_ZW_Chr3L':12, 'ZH_RAL_ZI_Chr3L':6, 'ZS_ZH_ZW_Chr2L':16, 'ZS_RAL_ZI_FR_SAfr_Chr2L':15483, 'ZW_RAL_ZI_Chr2L':16188, 'ZW_RAL_ZI_ChrX':15761, 'Zim_RAL_ZI_ChrX':15803, 'ZS_RAL_ZI_Chr3L':14940, 'ZS_RAL_ZI_autosome':58709, 'Zim_RAL_ZI_Chr2L':16120, 'ZW_RAL_ZI_autosome':29340, 'Zim_RAL_ZI_Chr3L':7, 'ZS_RAL_ZI_FR_SAfr_ChrX':15308, 'ZS_RAL_ZI_FR_SAfr_autosome':58175, 'Zim_RAL_ZI_Chr2R':13154, 'Zim_RAL_ZI_autosome':29325, 'ZS_RAL_ZI_ChrX':15322, 'ZS_ZH_ZW_autosome':11729, 'Zim_RAL_ZI_Chr3R':46, 'ZS_ZH_ZW_ChrX':15589, 'ZS_RAL_ZI_Chr2R':12513, 'ZW_RAL_ZI_Chr3R':8, 'ZH_RAL_ZI_ChrX':15319, 'ZS_RAL_ZI_FR_SAfr_Chr3R':15539, 'ZS_ZH_ZW_Chr3R':8, 'ZH_RAL_ZI_Chr3R':7, 'ZH_RAL_ZI_Chr2R':12199, 'ZS_ZH_ZW_Chr2R':11694, 'ZS_RAL_ZI_FR_SAfr_Chr2R':12252, 'ZW_RAL_ZI_Chr2R':13140, 'ZS_RAL_ZI_Chr3R':15755}
gene_sample_size = {'ZS_RAL_ZI_Chr2L':1514, 'ZH_RAL_ZI_autosome':1361, 'ZW_RAL_ZI_Chr3L':0, 'ZS_RAL_ZI_FR_SAfr_Chr3L':1363, 'ZS_ZH_ZW_Chr3L':1, 'ZH_RAL_ZI_Chr3L':1, 'ZS_ZH_ZW_Chr2L':3, 'ZS_RAL_ZI_FR_SAfr_Chr2L':1487, 'ZW_RAL_ZI_Chr2L':1478, 'ZW_RAL_ZI_ChrX':1180, 'Zim_RAL_ZI_ChrX':1183, 'ZS_RAL_ZI_Chr3L':1372, 'ZS_RAL_ZI_autosome':6114, 'Zim_RAL_ZI_Chr2L':1492, 'ZW_RAL_ZI_autosome':2861, 'Zim_RAL_ZI_Chr3L':0, 'ZS_RAL_ZI_FR_SAfr_ChrX':1214, 'ZS_RAL_ZI_FR_SAfr_autosome':6050, 'Zim_RAL_ZI_Chr2R':1469, 'Zim_RAL_ZI_autosome':2897, 'ZS_RAL_ZI_ChrX':1203, 'ZS_ZH_ZW_autosome':1423, 'Zim_RAL_ZI_Chr3R':8, 'ZS_ZH_ZW_ChrX':1253, 'ZS_RAL_ZI_Chr2R':1505, 'ZW_RAL_ZI_Chr3R':4, 'ZH_RAL_ZI_ChrX':1183, 'ZS_RAL_ZI_FR_SAfr_Chr3R':1725, 'ZS_ZH_ZW_Chr3R':5, 'ZH_RAL_ZI_Chr3R':5, 'ZH_RAL_ZI_Chr2R':1360, 'ZS_ZH_ZW_Chr2R':1423, 'ZS_RAL_ZI_FR_SAfr_Chr2R':1475, 'ZW_RAL_ZI_Chr2R':1458, 'ZS_RAL_ZI_Chr3R':1721}


#permutation function
def null_permutator(infile):

    print("reading dataframe for ", infile)

    #openning df
    df = pd.read_csv(infile)

    #filtering df
    genes_temp_df = df[df.FBgn != "None"]
    genes_df = genes_temp_df.drop_duplicates(subset = ["FBgn"])

    #setting counting lists
    gene_sites_count = []
    genes_represented_count = []
    FB_neurogenesis_sites_count = []
    FM_neurogenesis_sites_count = []
    FB_mating_behavior_sites_count = []
    FM_mating_behavior_sites_count = []
    FB_male_mating_behavior_sites_count = []
    FM_male_mating_behavior_sites_count = []
    FB_female_mating_behavior_sites_count = []
    FM_female_mating_behavior_sites_count = []

    FB_neurogenesis_genes_count = []
    FM_neurogenesis_genes_count = []
    FB_mating_behavior_genes_count = []
    FM_mating_behavior_genes_count = []
    FB_male_mating_behavior_genes_count = []
    FM_male_mating_behavior_genes_count = []
    FB_female_mating_behavior_genes_count = []
    FM_female_mating_behavior_genes_count = []

    print("finding sample sizes for ", infile)

    #loop to find sample sizes
    for i in site_sample_size:
        if i in infile:
            file_label = i
            site_sample = site_sample_size[i]
            gene_sample = gene_sample_size[i]
            break
        else:
            continue

    print("permutating...")

    #permutation
    for i in range(1000):
        
        #random sampling
        sampled_df = df.sample(n=site_sample)

        sampled_genes_df = genes_df.sample(n=gene_sample)

        #creating lists
        downsample_gene_sites = sampled_df['FBgn'].tolist()
        downsample_genes_represented = list(set(downsample_gene_sites))
        downsample_FB_neurogenesis_sites = sampled_df['FB_GO:0022008'].tolist()
        downsample_FM_neurogenesis_sites = sampled_df['FM_GO:0022008'].tolist()
        downsample_FB_mating_behavior_sites = sampled_df['FB_GO:0007617'].tolist()
        downsample_FM_mating_behavior_sites = sampled_df['FM_GO:0007617'].tolist()
        downsample_FB_male_mating_behavior_sites = sampled_df['FB_GO:0060179'].tolist()
        downsample_FM_male_mating_behavior_sites = sampled_df['FM_GO:0060179'].tolist()
        downsample_FB_female_mating_behavior_sites = sampled_df['FB_GO:0060180'].tolist()
        downsample_FM_female_mating_behavior_sites = sampled_df['FM_GO:0060180'].tolist()

        downsample_FB_neurogenesis_genes = sampled_genes_df['FB_GO:0022008'].tolist()
        downsample_FM_neurogenesis_genes = sampled_genes_df['FM_GO:0022008'].tolist()
        downsample_FB_mating_behavior_genes = sampled_genes_df['FB_GO:0007617'].tolist()
        downsample_FM_mating_behavior_genes = sampled_genes_df['FM_GO:0007617'].tolist()
        downsample_FB_male_mating_behavior_genes = sampled_genes_df['FB_GO:0060179'].tolist()
        downsample_FM_male_mating_behavior_genes = sampled_genes_df['FM_GO:0060179'].tolist()
        downsample_FB_female_mating_behavior_genes = sampled_genes_df['FB_GO:0060180'].tolist()
        downsample_FM_female_mating_behavior_genes = sampled_genes_df['FM_GO:0060180'].tolist()

        #counting and appending
        ##site counts

        #gene sites count
        gene_sites_count.append(len(downsample_gene_sites) - downsample_gene_sites.count("None"))

        #genes represented count
        genes_represented_count.append(len(downsample_genes_represented) - downsample_genes_represented.count("None"))

        #FB neurogenesis sites count
        FB_neurogenesis_sites_count.append(downsample_FB_neurogenesis_sites.count("Neurogenesis"))
        
        #FM neurogenesis sites count
        FM_neurogenesis_sites_count.append(downsample_FM_neurogenesis_sites.count("Neurogenesis"))
        
        #FB mating behavior sites count
        FB_mating_behavior_sites_count.append(downsample_FB_mating_behavior_sites.count("Mating_behavior"))
        
        #FM mating behavior sites count
        FM_mating_behavior_sites_count.append(downsample_FM_mating_behavior_sites.count("Mating_behavior"))
        
        #FB male mating behavior sites count
        FB_male_mating_behavior_sites_count.append(downsample_FB_male_mating_behavior_sites.count("Male_mating_behavior"))
        
        #FM male mating behavior sites count
        FM_male_mating_behavior_sites_count.append(downsample_FM_male_mating_behavior_sites.count("Male_mating_behavior"))
        
        #FB female mating behavior sites count
        FB_female_mating_behavior_sites_count.append(downsample_FB_female_mating_behavior_sites.count("Female_mating_behavior"))
        
        #FM female mating behavior sites count
        FM_female_mating_behavior_sites_count.append(downsample_FM_female_mating_behavior_sites.count("Female_mating_behavior"))

        ##gene counts
   
        #FB neurogenesis genes count
        FB_neurogenesis_genes_count.append(downsample_FB_neurogenesis_genes.count("Neurogenesis"))
        
        #FM neurogenesis genes count
        FM_neurogenesis_genes_count.append(downsample_FM_neurogenesis_genes.count("Neurogenesis"))
        
        #FB mating behavior genes count
        FB_mating_behavior_genes_count.append(downsample_FB_mating_behavior_genes.count("Mating_behavior"))
        
        #FM mating behavior genes count
        FM_mating_behavior_genes_count.append(downsample_FM_mating_behavior_genes.count("Mating_behavior"))
        
        #FB male mating behavior genes count
        FB_male_mating_behavior_genes_count.append(downsample_FB_male_mating_behavior_genes.count("Male_mating_behavior"))
        
        #FM male mating behavior genes count
        FM_male_mating_behavior_genes_count.append(downsample_FM_male_mating_behavior_genes.count("Male_mating_behavior"))
        
        #FB female mating behavior genes count
        FB_female_mating_behavior_genes_count.append(downsample_FB_female_mating_behavior_genes.count("Female_mating_behavior"))
        
        #FM female mating behavior genes count
        FM_female_mating_behavior_genes_count.append(downsample_FM_female_mating_behavior_genes.count("Female_mating_behavior"))



    #creating permutation dataframe
    permut_dict = {"num_of_gene_sites":gene_sites_count, "num_of_unique_genes":genes_represented_count, "Flybase_neurogenesis_sites":FB_neurogenesis_sites_count, "Flymine_neurogenesis_sites":FM_neurogenesis_sites_count, "Flybase_mating_behavior_sites":FB_mating_behavior_sites_count, "Flymine_mating_behavior_sites":FM_mating_behavior_sites_count, "Flybase_male_mating_behavior_sites":FB_male_mating_behavior_sites_count, "Flymine_male_mating_behavior_sites":FM_male_mating_behavior_sites_count, "Flybase_female_mating_behavior_sites":FB_female_mating_behavior_sites_count, "Flymine_female_mating_behavior_sites":FM_female_mating_behavior_sites_count, "Flybase_neurogenesis_genes":FB_neurogenesis_genes_count, "Flymine_neurogenesis_genes":FM_neurogenesis_genes_count, "Flybase_mating_behavior_genes":FB_mating_behavior_genes_count, "Flymine_mating_behavior_genes":FM_mating_behavior_genes_count, "Flybase_male_mating_behavior_genes":FB_male_mating_behavior_genes_count, "Flymine_male_mating_behavior_genes":FM_male_mating_behavior_genes_count, "Flybase_female_mating_behavior_genes":FB_female_mating_behavior_genes_count, "Flymine_female_mating_behavior_genes":FM_female_mating_behavior_genes_count}

    permut_df = pd.DataFrame.from_dict(permut_dict)

    #saving as csv
    permut_df.to_csv('null_permutation_{label}.csv'.format(label=file_label), index=True)

    print("new dataframe saved as: null_permutation_{label}.csv".format(label=file_label))






#for parallel mapping
csv_list = ['GO_sites_genes_null_Fst_ZS_RAL_ZI_Chr2L.csv', 'GO_sites_genes_null_Fst_ZH_RAL_ZI_autosome.csv', 'GO_sites_genes_null_Fst_ZW_RAL_ZI_Chr3L.csv', 'GO_sites_genes_null_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L.csv', 'GO_sites_genes_null_Fst_ZS_ZH_ZW_Chr3L.csv', 'GO_sites_genes_null_Fst_ZH_RAL_ZI_Chr3L.csv', 'GO_sites_genes_null_Fst_ZS_ZH_ZW_Chr2L.csv', 'GO_sites_genes_null_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L.csv', 'GO_sites_genes_null_Fst_ZW_RAL_ZI_Chr2L.csv', 'GO_sites_genes_null_Fst_ZW_RAL_ZI_ChrX.csv', 'GO_sites_genes_null_Fst_Zim_RAL_ZI_ChrX.csv', 'GO_sites_genes_null_Fst_ZS_RAL_ZI_Chr3L.csv', 'GO_sites_genes_null_Fst_ZS_RAL_ZI_autosome.csv', 'GO_sites_genes_null_Fst_Zim_RAL_ZI_Chr2L.csv', 'GO_sites_genes_null_Fst_ZW_RAL_ZI_autosome.csv', 'GO_sites_genes_null_Fst_Zim_RAL_ZI_Chr3L.csv', 'GO_sites_genes_null_Fst_ZS_RAL_ZI_FR_SAfr_ChrX.csv', 'GO_sites_genes_null_Fst_ZS_RAL_ZI_FR_SAfr_autosome.csv', 'GO_sites_genes_null_Fst_Zim_RAL_ZI_Chr2R.csv', 'GO_sites_genes_null_Fst_Zim_RAL_ZI_autosome.csv', 'GO_sites_genes_null_Fst_ZS_RAL_ZI_ChrX.csv', 'GO_sites_genes_null_Fst_ZS_ZH_ZW_autosome.csv', 'GO_sites_genes_null_Fst_Zim_RAL_ZI_Chr3R.csv', 'GO_sites_genes_null_Fst_ZS_ZH_ZW_ChrX.csv', 'GO_sites_genes_null_Fst_ZS_RAL_ZI_Chr2R.csv', 'GO_sites_genes_null_Fst_ZW_RAL_ZI_Chr3R.csv', 'GO_sites_genes_null_Fst_ZH_RAL_ZI_ChrX.csv', 'GO_sites_genes_null_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R.csv', 'GO_sites_genes_null_Fst_ZS_ZH_ZW_Chr3R.csv', 'GO_sites_genes_null_Fst_ZH_RAL_ZI_Chr3R.csv', 'GO_sites_genes_null_Fst_ZH_RAL_ZI_Chr2R.csv', 'GO_sites_genes_null_Fst_ZS_ZH_ZW_Chr2R.csv', 'GO_sites_genes_null_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R.csv', 'GO_sites_genes_null_Fst_ZW_RAL_ZI_Chr2R.csv', 'GO_sites_genes_null_Fst_ZS_RAL_ZI_Chr3R.csv']

#run in parallel
def run_in_parallel():
    pool = Pool(processes=8)
    pool.map(null_permutator, csv_list)


if __name__ == '__main__':
    run_in_parallel()


