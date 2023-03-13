import os, sys
import pandas as pd
import numpy as np

##gene name dictionary
#opening chrom maps
ChrX_map = pd.read_csv('../Dmel_genemap_ChrX.csv')
Chr2L_map = pd.read_csv('../Dmel_genemap_Chr2L.csv')
Chr2R_map = pd.read_csv('../Dmel_genemap_Chr2R.csv')
Chr3L_map = pd.read_csv('../Dmel_genemap_Chr3L.csv')
Chr3R_map = pd.read_csv('../Dmel_genemap_Chr3R.csv')

#making lists
fbgn_ls = ChrX_map['FBgn'].tolist() + Chr2L_map['FBgn'].tolist() + Chr2R_map['FBgn'].tolist() + Chr3L_map['FBgn'].tolist() + Chr3R_map['FBgn'].tolist()
gene_ls = ChrX_map['symbol'].tolist() + Chr2L_map['symbol'].tolist() + Chr2R_map['symbol'].tolist() + Chr3L_map['symbol'].tolist() + Chr3R_map['symbol'].tolist()
#making dictionary
gene_dict = dict(zip(fbgn_ls, gene_ls))
gene_dict.update({'':''})


##opening neuro list
#GO lists
GO_file = open('../neurogenesis_genes.txt', "r")
GO_genes = GO_file.read()
#to list
GO_list = GO_genes.split("\n")
GO_file.close()


#function
def get_gene_names_GO(subsetfile):
    """Function takes in a gene subset file from upset4genes.py and matches the FBgns to gene symbols and whether or not they are neurogenetic"""

    #openning subset file and making FBgns into list
    subset_file = open(subsetfile, "r")
    subsetstr = subset_file.read()
    subset = subsetstr.split("\n")  ###
    subset_file.close()

    #setting lists
    subsetsymbol = []
    neuro = []

    #looping through FBgns
    for i in subset:
        subsetsymbol.append(gene_dict[i])

        if i in GO_list and i != '':
            neuro.append("Neurogenesis")
        else:
            neuro.append("")

    #making dataframe
    dictionary = {"FBgn": subset, "Symbol": subsetsymbol, "Neurogenesis": neuro}
    df = pd.DataFrame(dictionary)

    #save
    df.to_csv(subsetfile.replace(".txt", "") + "_withsymbol.csv", index=False)




#running for all subset files
directory = '/Users/philipbaldassari/Desktop/annotated_persite_fst/top_percentage_csvs/upset_persite'

for file in os.listdir(directory):
    if file.endswith('.txt'):
        get_gene_names_GO(file)