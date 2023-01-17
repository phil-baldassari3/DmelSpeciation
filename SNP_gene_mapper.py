#script taking in list files of flybase point site coordinates and mapps them to genes

import os,sys
import pandas as pd
import numpy as np
from multiprocessing import Pool

#setting directory
directory = "/Users/philipbaldassari/Desktop/zim-cos_gene_GO_ann"
os.chdir(directory)


#get list of txt list files for parallel pool mapping
file_list = []

for file in os.listdir(directory):
    if file.endswith('.txt') and 'r6' in str(file):
        file_list.append(file)
    else:
        continue

print("\n########################################################################\n")
print(file_list)
print("\n########################################################################\n")


#gene mapper function
def gene_mapper(infile):


    #reading in Dmel gene maps
    ChrX_map = pd.read_csv('Dmel_genemap_ChrX.csv')
    Chr2L_map = pd.read_csv('Dmel_genemap_Chr2L.csv')
    Chr2R_map = pd.read_csv('Dmel_genemap_Chr2R.csv')
    Chr3L_map = pd.read_csv('Dmel_genemap_Chr3L.csv')
    Chr3R_map = pd.read_csv('Dmel_genemap_Chr3R.csv')

    print("Gene maps read. Starting gene mapping")

    print("starting the gene mapping on " + infile)



    #opening the coordinate file
    coordinate_file = open(infile, "r")
  
    #reading the file
    coordinate_data = coordinate_file.read()
  
    #coordinates to list
    coordinates = coordinate_data.split("\n")

    coordinate_file.close()
    print("coordinates from " + infile + " read into list")

    #gene lists
    FBgn_list = []
    gene_symbol_list = []

    #looping through coordinates
    for i in coordinates:

        if 'X' in i:
            df = ChrX_map
        elif '2L' in i:
            df = Chr2L_map
        elif '2R' in i:
            df = Chr2R_map
        elif '3L' in i:
            df = Chr3L_map
        elif '3R' in i:
            df = Chr3R_map
        else:
            print("Chromosome not found.")
            break
        
        #setting temp lists
        FBgn_list_temp = []
        gene_symbol_temp = []


        #looping through gene maps
        for row in df.itertuples():


            #tuple to list
            row_list = list(row)

            if int(i.split(':')[1]) >= row_list[4] and int(i.split(':')[1]) <= row_list[5]:
                FBgn_list_temp.append(row_list[1])
                gene_symbol_temp.append(str(row_list[2]))
            else:
                continue
        
        #temp lists to strings
        FBgn_str = "; ".join(FBgn_list_temp)
        gene_str = "; ".join(gene_symbol_temp)

        #checking if there is no hit
        if len(FBgn_str) == 0:
            FBgn_list.append('None')
            gene_symbol_list.append('None')
        else:
            FBgn_list.append(FBgn_str)
            gene_symbol_list.append(gene_str)
            print(gene_str)






    #making dataframe
    dictionary = {'FBgn' : FBgn_list, 'gene_symbol' : gene_symbol_list}

    gene_hits_df = pd.DataFrame(dictionary)

    #saving dataframe as csv
    gene_hits_df.to_csv('gene_hits_{file}.csv'.format(file=infile), index=False)

    print('csv of gene hits made from {file}'.format(file=infile))




#for parallele mapping

#run in parallel
def run_in_parallel():
    pool = Pool(processes=18)
    pool.map(gene_mapper, file_list)


if __name__ == '__main__':
    run_in_parallel()

