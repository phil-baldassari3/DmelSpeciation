#importing modules
import os, sys
import pandas as pd
import numpy as np
from statistics import mean
from statistics import stdev
import scipy.stats
import random
import matplotlib.pyplot as plt


#functions
def prepare4comparison(filename, GO_list_file):
    """Function takes in the islands csv filename and reads the correct gene map and neurogenesis list file"""

    #getting correct file
    if "ChrX" in filename:
        genemap_file = "gene_data/Dmel_genemap_ChrX.csv"
        chrom_len = 23542271
    elif "Chr2L" in filename:
        genemap_file = "gene_data/Dmel_genemap_Chr2L.csv"
        chrom_len = 23513712
    elif "Chr2R" in filename:
        genemap_file = "gene_data/Dmel_genemap_Chr2R.csv"
        chrom_len = 25286936
    elif "Chr3L" in filename:
        genemap_file = "gene_data/Dmel_genemap_Chr3L.csv"
        chrom_len = 28110227
    elif "Chr3R" in filename:
        genemap_file = "gene_data/Dmel_genemap_Chr3R.csv"
        chrom_len = 32079331
    else:
        print("SOMETHING IS WRONG")


    #gene map to df
    genemap_df = pd.read_csv(genemap_file)

    #enrichment list file
    with open("gene_data/" + GO_list_file, 'r') as GOfile:
        GOstr = GOfile.read()
        GO_list = GOstr.split('\n')

    return chrom_len, genemap_df, GO_list


def count_from_islands(comparison_name, islands_list, genemap, GOlist):
    """Fucntion uses the island_list and counts the number of GO term genes in the island"""

    #lists to be filled
    GO_counts = []
    islands_len = []

    islands = [eval(x) for x in islands_list if type(x) != float]

    #looping through islands
    for i in islands:
        
        #filtering gene map for genes containing in the island
        genemap_filtered = genemap[
            (genemap['start'] <= i[0]) & (genemap['stop'] >= i[0]) | 
            (genemap['start'] >= i[0]) & (genemap['stop'] <= i[1]) | 
            (genemap['start'] <= i[1]) & (genemap['stop'] >= i[1])
            ]
        
        genes = genemap_filtered["FBgn"].to_list()

        #counting how many GO term genes
        count = 0
        for g in genes:
            if g in GOlist:
                count += 1
            else:
                continue

        #appending to lists
        GO_counts.append(count)
        islands_len.append(i[1] - i[0])


    return GO_counts, islands_len


def enrichment_per_island(GOcounts, is_lengths, gene_map, chromosome_length):
    """Function uses the GO counts, island lengths, gene map, and chromosome length to perform a Z-test for enruichment of the GO term in the island"""