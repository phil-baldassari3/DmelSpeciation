#script takes the top percentage subset of the loci table and outputs a dataframe of the genes contained in the set with SNP and Fst information

#importing modules
import os,sys
import pandas as pd
import numpy as np
from multiprocessing import Pool


#functions

def gene_ls_set(ls):
    """Function takes a list of genes (use .to_list() to extract df column) in the form [gene1, gene2; gene1, ...] and outputs a cleaned list of genes and a list of unique gene names from the full list"""

    ls = [str(i) for i in ls]

    strls = "; ".join(ls)
    newls = strls.split("; ")

    setls = list(set(newls))

    return setls

   


def annotate_GO(dataframe, GO_term, GO_listfile):
    """
    This method takes three parameters:
    dataframe: df of genes
    GO_term: string of the GO term you whan to annotate.
    GO_listfile: string of a filename that contains a list of genes, one gene per line, that fir your chosen GO term
    """
    
    ##GO lists
    GO_file = open(GO_listfile, "r")
    GO_genes = GO_file.read()
    #to list
    GO_list = GO_genes.split("\n")
    GO_file.close()

    
    #gene list from table
    gene_list = dataframe['FBgn'].tolist()

    

    #setting list
    go = []

    #looping through null genes
    for gene in gene_list:

        if gene in GO_list:
            go.append(GO_term)
        else:
            go.append("None")
        

    #annotating table
    dictionary = {GO_term : go}
    GO_df = pd.DataFrame(dictionary)

    newdf = pd.merge(dataframe, GO_df, left_index=True, right_index=True)

    return newdf


def gene_df_maker(locitable):
    """Loops through genes frojm the loci table and computes the number of SNPs per gene and average Fst of that gene. Returns a new df. This function used gene_ls_set()."""

    #opening chrom maps
    ChrX_map = pd.read_csv('Dmel_genemap_ChrX.csv')
    Chr2L_map = pd.read_csv('Dmel_genemap_Chr2L.csv')
    Chr2R_map = pd.read_csv('Dmel_genemap_Chr2R.csv')
    Chr3L_map = pd.read_csv('Dmel_genemap_Chr3L.csv')
    Chr3R_map = pd.read_csv('Dmel_genemap_Chr3R.csv')

    #making lists
    fbgn_ls = ChrX_map['FBgn'].tolist() + Chr2L_map['FBgn'].tolist() + Chr2R_map['FBgn'].tolist() + Chr3L_map['FBgn'].tolist() + Chr3R_map['FBgn'].tolist()
    gene_ls = ChrX_map['symbol'].tolist() + Chr2L_map['symbol'].tolist() + Chr2R_map['symbol'].tolist() + Chr3L_map['symbol'].tolist() + Chr3R_map['symbol'].tolist()
    #making dictionary
    gene_dict = dict(zip(fbgn_ls, gene_ls))




    #gene columns to lists
    raw_fbgnls = locitable["FBgn"].to_list()

    #getting a list of unique genes
    setfbgns = gene_ls_set(raw_fbgnls)
    setfbgns.remove("None")

    #lists for new df
    fbgn = []
    gene = []
    num = []
    avg = []

    #looping through genes and filtering dataframe for computation
    for i in setfbgns:

        genedf = locitable[locitable['FBgn'].str.contains(i)]

        #appending to lists
        fbgn.append(i)
        gene.append(gene_dict[i])
        num.append(len(genedf))
        avg.append(genedf['Avg_Fst'].mean())

    #new df
    dictionary = {"FBgn": fbgn, "gene_symbol": gene, "num_of_SNPs": num, "Avg_Fst": avg}
    newdf = pd.DataFrame.from_dict(dictionary)

    return newdf





def main_func(infile):
    """main function of this program. Takes an infile (csv) and outputs a new gene df adn saves it as a csv."""

    df = pd.read_csv(infile)

    #new gene df
    genedf = gene_df_maker(df)

    #annotating gene df
    finalgenedf = annotate_GO(genedf, "Neurogenesis", "neurogenesis_genes.txt")

    #df to file
    finalgenedf.to_csv(infile.replace('annotated', 'genes'), index=False)











directory = "/Users/philipbaldassari/Desktop/annotated_persite_fst/top_percentage_csvs"
os.chdir(directory)

#get list of txt list files for parallel pool mapping
file_list = []

for file in os.listdir(directory):
    if file.endswith('.csv') and 'annotated' in file:
        file_list.append(file)
    else:
        continue

print(len(file_list))




#run in parallel
def run_in_parallel():
    pool = Pool(processes=18)
    pool.map(main_func, file_list)


if __name__ == '__main__':
    run_in_parallel()





