#script to annotate the null genes for neurogenesis and mating behavior GO terms and child terms
import os,sys
import pandas as pd
import numpy as np
from multiprocessing import Pool

#setting directory
directory = "/Users/philipbaldassari/Desktop/zim-cos_downsampled/null_condition/top0.01_gene_sites_genes"
os.chdir(directory)

'''
#get list of txt list files for parallel pool mapping
file_list = []

for file in os.listdir(directory):
    if file.endswith('.txt') and 'GO' not in str(file):
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


#GO counter
def GO_counter(infile):

    gene_sites_file = open(infile, "r")
    gene_sites_data = gene_sites_file.read()
    #to lists
    gene_sites_list = gene_sites_data.split("\n")

    gene_sites_file.close()

    unique_genes_list = list(set(gene_sites_list))

    print("opened ", infile, " and saved data to lists")

    #setting lists
    GO_neurogenesis_sites_FB = []
    GO_mating_behavior_sites_FB = []
    GO_male_mating_behavior_sites_FB = []
    GO_female_mating_behavior_sites_FB = []
    GO_neurogenesis_sites_FM = []
    GO_mating_behavior_sites_FM = []
    GO_male_mating_behavior_sites_FM = []
    GO_female_mating_behavior_sites_FM = []


    print("looping through gene sites of ", infile)

    #looping through top 1% gene sites
    for gene in gene_sites_list:
        
        if gene in FB_neurogenesis_GO:
            GO_neurogenesis_sites_FB.append('Neurogenesis')
        else:
            GO_neurogenesis_sites_FB.append('None')
        

        if gene in FB_mating_behavior_GO:
            GO_mating_behavior_sites_FB.append('Mating_behavior')
        else:
            GO_mating_behavior_sites_FB.append('None')
        

        if gene in FB_male_mating_behavior_GO:
            GO_male_mating_behavior_sites_FB.append('Male_mating_behavior')
        else:
            GO_male_mating_behavior_sites_FB.append('None')
        

        if gene in FB_female_mating_behavior_GO:
            GO_female_mating_behavior_sites_FB.append('Female_mating_behavior')
        else:
            GO_female_mating_behavior_sites_FB.append('None')
        


        if gene in FM_neurogenesis_GO:
            GO_neurogenesis_sites_FM.append('Neurogenesis')
        else:
            GO_neurogenesis_sites_FM.append('None')
        

        if gene in FM_mating_behavior_GO:
            GO_mating_behavior_sites_FM.append('Mating_behavior')
        else:
            GO_mating_behavior_sites_FM.append('None')
        

        if gene in FM_male_mating_behavior_GO:
            GO_male_mating_behavior_sites_FM.append('Male_mating_behavior')
        else:
            GO_male_mating_behavior_sites_FM.append('None')
        

        if gene in FM_female_mating_behavior_GO:
            GO_female_mating_behavior_sites_FM.append('Female_mating_behavior')
        else:
            GO_female_mating_behavior_sites_FM.append('None')


    #counting number of gene sites
    num_gene_sites = len(gene_sites_list)
    num_unique_genes = len(unique_genes_list)

    num_neurogenesis_sites_FB = GO_neurogenesis_sites_FB.count("Neurogenesis")
    num_neurogenesis_sites_FM = GO_neurogenesis_sites_FM.count("Neurogenesis")
    num_mating_behavior_sites_FB = GO_mating_behavior_sites_FB.count("Mating_behavior")
    num_mating_behavior_sites_FM = GO_mating_behavior_sites_FM.count("Mating_behavior")
    num_male_mating_behavior_sites_FB = GO_male_mating_behavior_sites_FB.count("Male_mating_behavior")
    num_male_mating_behavior_sites_FM = GO_male_mating_behavior_sites_FM.count("Male_mating_behavior")
    num_female_mating_behavior_sites_FB = GO_female_mating_behavior_sites_FB.count("Female_mating_behavior")
    num_female_mating_behavior_sites_FM = GO_female_mating_behavior_sites_FM.count("Female_mating_behavior")

    print("done counting of gene sites for ", infile)

    #setting lists
    GO_neurogenesis_genes_FB = []
    GO_mating_behavior_genes_FB = []
    GO_male_mating_behavior_genes_FB = []
    GO_female_mating_behavior_genes_FB = []
    GO_neurogenesis_genes_FM = []
    GO_mating_behavior_genes_FM = []
    GO_male_mating_behavior_genes_FM = []
    GO_female_mating_behavior_genes_FM = []

    print('looping though unique genes of ', infile)

    #looping through top 1% genes
    for gene in unique_genes_list:
        
        if gene in FB_neurogenesis_GO:
            GO_neurogenesis_genes_FB.append('Neurogenesis')
        else:
            GO_neurogenesis_genes_FB.append('None')
        

        if gene in FB_mating_behavior_GO:
            GO_mating_behavior_genes_FB.append('Mating_behavior')
        else:
            GO_mating_behavior_genes_FB.append('None')
        

        if gene in FB_male_mating_behavior_GO:
            GO_male_mating_behavior_genes_FB.append('Male_mating_behavior')
        else:
            GO_male_mating_behavior_genes_FB.append('None')
        

        if gene in FB_female_mating_behavior_GO:
            GO_female_mating_behavior_genes_FB.append('Female_mating_behavior')
        else:
            GO_female_mating_behavior_genes_FB.append('None')
        



        if gene in FM_neurogenesis_GO:
            GO_neurogenesis_genes_FM.append('Neurogenesis')
        else:
            GO_neurogenesis_genes_FM.append('None')
        

        if gene in FM_mating_behavior_GO:
            GO_mating_behavior_genes_FM.append('Mating_behavior')
        else:
            GO_mating_behavior_genes_FM.append('None')
        

        if gene in FM_male_mating_behavior_GO:
            GO_male_mating_behavior_genes_FM.append('Male_mating_behavior')
        else:
            GO_male_mating_behavior_genes_FM.append('None')
        

        if gene in FM_female_mating_behavior_GO:
            GO_female_mating_behavior_genes_FM.append('Female_mating_behavior')
        else:
            GO_female_mating_behavior_genes_FM.append('None')



    num_neurogenesis_genes_FB = GO_neurogenesis_genes_FB.count("Neurogenesis")
    num_neurogenesis_genes_FM = GO_neurogenesis_genes_FM.count("Neurogenesis")
    num_mating_behavior_genes_FB = GO_mating_behavior_genes_FB.count("Mating_behavior")
    num_mating_behavior_genes_FM = GO_mating_behavior_genes_FM.count("Mating_behavior")
    num_male_mating_behavior_genes_FB = GO_male_mating_behavior_genes_FB.count("Male_mating_behavior")
    num_male_mating_behavior_genes_FM = GO_male_mating_behavior_genes_FM.count("Male_mating_behavior")
    num_female_mating_behavior_genes_FB = GO_female_mating_behavior_genes_FB.count("Female_mating_behavior")
    num_female_mating_behavior_genes_FM = GO_female_mating_behavior_genes_FM.count("Female_mating_behavior")

    print('done counts of unique genes for ', infile)
    

    with open("count_{file}".format(file=infile), "w") as counts_file:
        counts_file.write(str(infile) + "\n"+ "\n")
        counts_file.write("num_gene_sites: " + str(num_gene_sites) + "\n")
        counts_file.write("num_unique_genes: " + str(num_unique_genes) + "\n")
        counts_file.write("num_neurogenesis_sites_FB: " + str(num_neurogenesis_sites_FB) + "\n")
        counts_file.write("num_neurogenesis_sites_FM: " + str(num_neurogenesis_sites_FM) + "\n")
        counts_file.write("num_mating_behavior_sites_FB: " + str(num_mating_behavior_sites_FB) + "\n")
        counts_file.write("num_mating_behavior_sites_FM: " + str(num_mating_behavior_sites_FM) + "\n")
        counts_file.write("num_male_mating_behavior_sites_FB: " + str(num_male_mating_behavior_sites_FB) + "\n")
        counts_file.write("num_male_mating_behavior_sites_FM: " + str(num_male_mating_behavior_sites_FM) + "\n")
        counts_file.write("num_female_mating_behavior_sites_FB: " + str(num_female_mating_behavior_sites_FB) + "\n")
        counts_file.write("num_female_mating_behavior_sites_FM: " + str(num_female_mating_behavior_sites_FM) + "\n" + "\n")
        counts_file.write("num_neurogenesis_genes_FB: " + str(num_neurogenesis_genes_FB) + "\n")
        counts_file.write("num_neurogenesis_genes_FM: " + str(num_neurogenesis_genes_FM) + "\n")
        counts_file.write("num_mating_behavior_genes_FB: " + str(num_mating_behavior_genes_FB) + "\n")
        counts_file.write("num_mating_behavior_genes_FM: " + str(num_mating_behavior_genes_FM) + "\n")
        counts_file.write("num_male_mating_behavior_genes_FB: " + str(num_male_mating_behavior_genes_FB) + "\n")
        counts_file.write("num_male_mating_behavior_genes_FM: " + str(num_male_mating_behavior_genes_FM) + "\n")
        counts_file.write("num_female_mating_behavior_genes_FB: " + str(num_female_mating_behavior_genes_FB) + "\n")
        counts_file.write("num_female_mating_behavior_genes_FM: " + str(num_female_mating_behavior_genes_FM))


    print('data for ', infile, 'writting to .txt file')

#for parallel mapping
csv_list = ['top0.01_gene_sites_ZH_RAL_ZI_Chr2R.txt', 'top0.01_gene_sites_ZS_ZH_ZW_Chr3L.txt', 'top0.01_gene_sites_ZW_RAL_ZI_Chr2R.txt', 'top0.01_gene_sites_Zim_RAL_ZI_Chr3L.txt', 'top0.01_gene_sites_ZS_RAL_ZI_Chr3R.txt', 'top0.01_gene_sites_ZS_RAL_ZI_FR_SAfr_Chr3R.txt', 'top0.01_gene_sites_ZS_RAL_ZI_FR_SAfr_Chr2R.txt', 'top0.01_gene_sites_ZS_RAL_ZI_Chr2R.txt', 'top0.01_gene_sites_Zim_RAL_ZI_Chr2L.txt', 'top0.01_gene_sites_ZW_RAL_ZI_ChrX.txt', 'top0.01_gene_sites_ZW_RAL_ZI_Chr3R.txt', 'top0.01_gene_sites_ZS_ZH_ZW_Chr2L.txt', 'top0.01_gene_sites_ZH_RAL_ZI_Chr3R.txt', 'top0.01_gene_sites_Zim_RAL_ZI_autosomes.txt', 'top0.01_gene_sites_ZS_RAL_ZI_FR_SAfr_autosomes.txt', 'top0.01_gene_sites_ZS_ZH_ZW_ChrX.txt', 'top0.01_gene_sites_ZH_RAL_ZI_autosomes.txt', 'top0.01_gene_sites_ZS_ZH_ZW_autosomes.txt', 'top0.01_gene_sites_ZS_RAL_ZI_autosomes.txt', 'top0.01_gene_sites_ZS_RAL_ZI_FR_SAfr_ChrX.txt', 'top0.01_gene_sites_ZS_RAL_ZI_ChrX.txt', 'top0.01_gene_sites_ZS_ZH_ZW_Chr3R.txt', 'top0.01_gene_sites_ZH_RAL_ZI_ChrX.txt', 'top0.01_gene_sites_ZW_RAL_ZI_Chr2L.txt', 'top0.01_gene_sites_Zim_RAL_ZI_Chr3R.txt', 'top0.01_gene_sites_ZS_RAL_ZI_Chr3L.txt', 'top0.01_gene_sites_Zim_RAL_ZI_ChrX.txt', 'top0.01_gene_sites_ZS_RAL_ZI_FR_SAfr_Chr3L.txt', 'top0.01_gene_sites_ZS_RAL_ZI_FR_SAfr_Chr2L.txt', 'top0.01_gene_sites_ZS_RAL_ZI_Chr2L.txt', 'top0.01_gene_sites_Zim_RAL_ZI_Chr2R.txt', 'top0.01_gene_sites_ZW_RAL_ZI_Chr3L.txt', 'top0.01_gene_sites_ZW_RAL_ZI_autosomes.txt', 'top0.01_gene_sites_ZS_ZH_ZW_Chr2R.txt', 'top0.01_gene_sites_ZH_RAL_ZI_Chr3L.txt']

#run in parallel
def run_in_parallel():
    pool = Pool(processes=8)
    pool.map(GO_counter, csv_list)


if __name__ == '__main__':
    run_in_parallel()


       

