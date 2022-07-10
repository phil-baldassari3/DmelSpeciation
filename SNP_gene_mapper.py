#script taking in list files of flybase point site coordinates and mapps them to genes

import os,sys
import pandas as pd
import numpy as np
from multiprocessing import Pool

#setting directory
directory = "set"
os.chdir(directory)

'''
#get list of txt list files for parallel pool mapping
file_list = []

for file in os.listdir(directory):
    if file.endswith('.txt'):
        file_list.append(file)
    else:
        continue

print(file_list)
'''

#reading in Dmel gene maps
ChrX_map = pd.read_csv('Dmel_genemap_ChrX.csv')
Chr2L_map = pd.read_csv('Dmel_genemap_Chr2L.csv')
Chr2R_map = pd.read_csv('Dmel_genemap_Chr2R.csv')
Chr3L_map = pd.read_csv('Dmel_genemap_Chr3L.csv')
Chr3R_map = pd.read_csv('Dmel_genemap_Chr3R.csv')

print("Gene maps read. Starting gene mapping")

wait = input("Press Enter to continue.")


#gene mapper function
def gene_mapper(infile):

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
            continue
        
        #setting counter and dummy list
        counter = 0
        dummy_list = []

        #looping through gene maps
        for row in df.itertuples():

            counter += 1

            #tuple to list
            row_list = list(row)

            if int(i.split(':')[1]) >= row_list[4] and int(i.split(':')[1]) <= row_list[5]:
                FBgn_list.append(row_list[1])
                gene_symbol_list.append(row_list[2])
                print("Hit!")
                break
            else:
                dummy_list.append(0)
            
        #checking if there is no hit
        if counter == len(dummy_list):
            FBgn_list.append('None')
            gene_symbol_list.append('None')
        else:
            continue

    #making dataframe
    dictionary = {'FBgn' : FBgn_list, 'gene_symbol' : gene_symbol_list}

    gene_hits_df = pd.DataFrame(dictionary)

    #saving dataframe as csv
    gene_hits_df.to_csv('gene_hits_{file}.csv'.format(file=infile), index=False)

    print('csv of gene hits made from {file}'.format(file=infile))




#for parallele mapping
txt_list = []

#run in parallel
def run_in_parallel():
    pool = Pool(processes=8)
    pool.map(gene_mapper, txt_list)


if __name__ == '__main__':
    run_in_parallel()

