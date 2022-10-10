#Script to permutate over the null set to count relevantly functional genes

#importing modules
import os,sys
import pandas as pd
import numpy as np
from multiprocessing import Pool
import random

#setting directory
directory = "/Users/philipbaldassari/Desktop/zim-cos_gene_GO_ann"
os.chdir(directory)

'''
#get list of txt list files for parallel pool mapping
file_list = []

for file in os.listdir(directory):
    if file.endswith('.csv') and 'GO_sites_genes' in str(file) and "top" not in str(file):
        file_list.append(file)
    else:
        continue

print(file_list)
'''



#permutation function
def null_permutator(infile):

    #1% sample sizes
    site_sample_size = {'FR_vs_RAL_ZI_SAfr_Fst_ChrX':7963, 'FR_vs_RAL_ZI_SAfr_Fst_autosomes':30484, 'RAL_vs_FR_ZI_SAfr_Fst_ChrX':7963, 'RAL_vs_FR_ZI_SAfr_Fst_autosomes':30484, 'SAfr_vs_RAL_FR_ZI_Fst_ChrX':7963, 'SAfr_vs_RAL_FR_ZI_Fst_autosomes':30484, 'ZH_RAL_ZI_Fst_ChrX':7660, 'ZH_RAL_ZI_Fst_autosomes':6106, 'ZI_vs_RAL_FR_SAfr_Fst_ChrX':7963, 'ZI_vs_RAL_FR_SAfr_Fst_autosomes':30484, 'ZS_RAL_ZI_FR_SAfr_Fst_ChrX':7654, 'ZS_RAL_ZI_FR_SAfr_Fst_autosomes':29088, 'ZS_RAL_ZI_Fst_ChrX':7661, 'ZS_RAL_ZI_Fst_autosomes':29355, 'ZS_ZH_ZW_Fst_ChrX':7795, 'ZS_ZH_ZW_Fst_autosomes':5865, 'ZW_RAL_ZI_Fst_ChrX':7881, 'ZW_RAL_ZI_Fst_autosomes':14670}


    ##dictionary for gene lengths
    #opening chrom maps
    ChrX_map = pd.read_csv('Dmel_genemap_ChrX.csv')
    Chr2L_map = pd.read_csv('Dmel_genemap_Chr2L.csv')
    Chr2R_map = pd.read_csv('Dmel_genemap_Chr2R.csv')
    Chr3L_map = pd.read_csv('Dmel_genemap_Chr3L.csv')
    Chr3R_map = pd.read_csv('Dmel_genemap_Chr3R.csv')
    #making lists
    gene_ls = ChrX_map['FBgn'].tolist() + Chr2L_map['FBgn'].tolist() + Chr2R_map['FBgn'].tolist() + Chr3L_map['FBgn'].tolist() + Chr3R_map['FBgn'].tolist()
    gene_len_ls =ChrX_map['length'].tolist() + Chr2L_map['length'].tolist() + Chr2R_map['length'].tolist() + Chr3L_map['length'].tolist() + Chr3R_map['length'].tolist()
    #making dictionaries
    genes_length = dict(zip(gene_ls, gene_len_ls))
    genes_length.update({'': 1})
    print("made gene length dictionary")



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




############################################################################################################################################




    print("reading dataframe for ", infile)

    #openning df
    df = pd.read_csv(infile)


    #setting counting lists
    

    print("finding sample sizes for ", infile)

    #loop to find sample sizes
    for i in site_sample_size:
        if i in infile:
            file_label = i
            site_sample = site_sample_size[i]
            print(i, site_sample)
            break
        else:
            continue

    print("permutating...")

    #setting lists
    neurogenesis =[]
    mating_behavior = []
    male_mating_behavior = []
    female_mating_behvaior = []

    weighted_neurogenesis =[]
    weighted_mating_behavior = []
    weighted_male_mating_behavior = []
    weighted_female_mating_behvaior = []


    #permutation
    for i in range(100000):
        
        #random sampling
        sampled_df = df.sample(n=site_sample)


        #Neurogenesis
        neurogenesis_df = sampled_df[sampled_df.Neurogenesis == "Neurogenesis"]

        neurogenesis_gene_temporary = neurogenesis_df['FBgn'].tolist()
        neurogenesis_gene_temp = [i for i in neurogenesis_gene_temporary if i != "None"]
        neurogenesis_gene_str = "; ".join(neurogenesis_gene_temp)
        full_neurogenesis_gene_ls = neurogenesis_gene_str.split("; ")  ###

        neurogenesis_unique_gene_ls_temp = list(set(full_neurogenesis_gene_ls))
        neurogenesis_unique_gene_ls = [j for j in neurogenesis_unique_gene_ls_temp if j in FB_neurogenesis_GO or j in FM_neurogenesis_GO]
        
        if '' in neurogenesis_unique_gene_ls:
            neurogenesis_unique_gene_ls.remove('') ###
        else:
            neurogenesis_unique_gene_ls = neurogenesis_unique_gene_ls

        #unweighted
        neurogenesis.append(len(neurogenesis_unique_gene_ls))

        #weighted
        neurogenesis_SNP_den = []

        for gene in neurogenesis_unique_gene_ls:
            gene_len = genes_length[gene]

            if len(neurogenesis_unique_gene_ls) == 0:
                neurogenesis_SNP_den.append(0)
            else:
                neurogenesis_SNP_den.append(((full_neurogenesis_gene_ls.count(gene)) / (gene_len)) * 1000)

        weighted_neurogenesis.append(sum(neurogenesis_SNP_den))

        #Mating Behavior
        mating_behavior_df = sampled_df[sampled_df.Mating_behavior == "Mating_behavior"]

        mating_behavior_gene_temporary = mating_behavior_df['FBgn'].tolist()
        mating_behavior_gene_temp = [i for i in mating_behavior_gene_temporary if i != "None"]
        mating_behavior_gene_str = "; ".join(mating_behavior_gene_temp)
        full_mating_behavior_gene_ls = mating_behavior_gene_str.split("; ")  ###

        mating_unique_gene_ls_temp = list(set(full_mating_behavior_gene_ls))
        mating_unique_gene_ls = [j for j in mating_unique_gene_ls_temp if j in FB_mating_behavior_GO or j in FM_mating_behavior_GO]
        
        if '' in mating_unique_gene_ls:
            mating_unique_gene_ls.remove('') ###
        else:
            mating_unique_gene_ls = mating_unique_gene_ls


        #unweighted
        mating_behavior.append(len(mating_unique_gene_ls))

        #weighted
        mating_behavior_SNP_den = []

        for gene in mating_unique_gene_ls:
            gene_len = genes_length[gene]

            if len(mating_unique_gene_ls) == 0:
                mating_behavior_SNP_den.append(0)
            else:
                mating_behavior_SNP_den.append(((full_mating_behavior_gene_ls.count(gene)) / (gene_len)) * 1000)

        weighted_mating_behavior.append(sum(mating_behavior_SNP_den))


        #Male mating behavior
        male_mating_behavior_df = sampled_df[sampled_df.Male_mating_behavior == "Male_mating_behavior"]

        male_mating_behavior_gene_temporary = male_mating_behavior_df['FBgn'].tolist()
        male_mating_behavior_gene_temp = [i for i in male_mating_behavior_gene_temporary if i != "None"]
        male_mating_behavior_gene_str = "; ".join(male_mating_behavior_gene_temp)
        full_male_mating_behavior_gene_ls = male_mating_behavior_gene_str.split("; ")  ###

        male_mating_unique_gene_ls_temp = list(set(full_male_mating_behavior_gene_ls))
        male_mating_unique_gene_ls = [j for j in male_mating_unique_gene_ls_temp if j in FB_male_mating_behavior_GO or j in FM_male_mating_behavior_GO]
        
        if '' in male_mating_unique_gene_ls:
            male_mating_unique_gene_ls.remove('') ###
        else:
            male_mating_unique_gene_ls = male_mating_unique_gene_ls

        #unweighted
        male_mating_behavior.append(len(male_mating_unique_gene_ls))

        #weighted
        male_mating_behavior_SNP_den = []

        for gene in male_mating_unique_gene_ls:
            gene_len = genes_length[gene]

            if len(male_mating_unique_gene_ls) == 0:
                male_mating_behavior_SNP_den.append(0)
            else:
                male_mating_behavior_SNP_den.append(((full_male_mating_behavior_gene_ls.count(gene)) / (gene_len)) * 1000)

        weighted_male_mating_behavior.append(sum(male_mating_behavior_SNP_den))


        #Female mating behavior
        female_mating_behavior_df = sampled_df[sampled_df.Female_mating_behavior == "female_mating_behavior"]

        female_mating_behavior_gene_temporary = female_mating_behavior_df['FBgn'].tolist()
        female_mating_behavior_gene_temp = [i for i in female_mating_behavior_gene_temporary if i != "None"]
        female_mating_behavior_gene_str = "; ".join(female_mating_behavior_gene_temp)
        full_female_mating_behavior_gene_ls = female_mating_behavior_gene_str.split("; ")  ###

        female_mating_unique_gene_ls_temp = list(set(full_female_mating_behavior_gene_ls))
        female_mating_unique_gene_ls = [j for j in female_mating_unique_gene_ls_temp if (j in FB_female_mating_behavior_GO) | (j in FM_female_mating_behavior_GO)]
        
        if '' in female_mating_unique_gene_ls:
            female_mating_unique_gene_ls.remove('') ###
        else:
            female_mating_unique_gene_ls = female_mating_unique_gene_ls

        #unweighted
        female_mating_behvaior.append(len(female_mating_unique_gene_ls))

        #weighted
        female_mating_behavior_SNP_den = []

        for gene in female_mating_unique_gene_ls:
            gene_len = genes_length[gene]

            if len(female_mating_unique_gene_ls) == 0:
                female_mating_behavior_SNP_den.append(0)
            else:
                female_mating_behavior_SNP_den.append(((full_female_mating_behavior_gene_ls.count(gene)) / (gene_len)) * 1000)

        weighted_female_mating_behvaior.append(sum(female_mating_behavior_SNP_den))

        


    #creating permutation dataframe
    permut_dict = {"Neurogenesis":neurogenesis, "Mating_behavior": mating_behavior, "Male_mating_behavior": male_mating_behavior, "Female_mating_behavior": female_mating_behvaior, "Weighted_Neurogenesis":weighted_neurogenesis, "Weighted_Mating_behavior": weighted_mating_behavior, "Weighted_Male_mating_behavior": weighted_male_mating_behavior, "Weighted_Female_mating_behavior": weighted_female_mating_behvaior}

    permut_df = pd.DataFrame.from_dict(permut_dict)

    #saving as csv
    permut_df.to_csv('null_permutation_0.5percent_{label}.csv'.format(label=file_label), index=True)

    print("new dataframe saved as: null_permutation_{label}.csv".format(label=file_label))





#for parallel mapping
csv_list = ['GO_sites_genes_SAfr_vs_RAL_FR_ZI_Fst_autosomes.csv', 'GO_sites_genes_ZS_ZH_ZW_Fst_ChrX.csv', 'GO_sites_genes_FR_vs_RAL_ZI_SAfr_Fst_autosomes.csv', 'GO_sites_genes_FR_vs_RAL_ZI_SAfr_Fst_ChrX.csv', 'GO_sites_genes_ZW_RAL_ZI_Fst_autosomes.csv', 'GO_sites_genes_ZS_RAL_ZI_FR_SAfr_Fst_ChrX.csv', 'GO_sites_genes_ZS_RAL_ZI_Fst_autosomes.csv', 'GO_sites_genes_ZS_RAL_ZI_FR_SAfr_Fst_autosomes.csv', 'GO_sites_genes_ZI_vs_RAL_FR_SAfr_Fst_autosomes.csv', 'GO_sites_genes_ZI_vs_RAL_FR_SAfr_Fst_ChrX.csv', 'GO_sites_genes_ZH_RAL_ZI_Fst_ChrX.csv', 'GO_sites_genes_ZW_RAL_ZI_Fst_ChrX.csv', 'GO_sites_genes_SAfr_vs_RAL_FR_ZI_Fst_ChrX.csv', 'GO_sites_genes_ZH_RAL_ZI_Fst_autosomes.csv', 'GO_sites_genes_ZS_ZH_ZW_Fst_autosomes.csv', 'GO_sites_genes_RAL_vs_FR_ZI_SAfr_Fst_autosomes.csv', 'GO_sites_genes_RAL_vs_FR_ZI_SAfr_Fst_ChrX.csv', 'GO_sites_genes_ZS_RAL_ZI_Fst_ChrX.csv']

#run in parallel
def run_in_parallel():
    pool = Pool(processes=18)
    pool.map(null_permutator, csv_list)


if __name__ == '__main__':
    run_in_parallel()



