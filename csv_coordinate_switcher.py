#script to replace the r5 coordnates with r6 coordinates in the csv files

#importing modules
import os,sys
import pandas as pd
import numpy as np
from multiprocessing import Pool

#setting directory
directory = "/Users/philipbaldassari/Desktop/zim-cos_downsampled/sliding_window_fst"
os.chdir(directory)

'''
#get list of txt list files for parallel pool mapping
file_list = []

for file in os.listdir(directory):
    if file.endswith('.csv'):
        file_list.append(file)
    else:
        continue

print(file_list)
'''

def coordinate_switcher(infile):

    #opening dataframe
    df = pd.read_csv(infile)

    print("opened ", infile)

    df.drop('Site', axis=1, inplace=True)

    print("dropped r5 sites from ", infile)

    for file in os.listdir(directory):
        if file.endswith('.txt') and infile.split('.')[0] in file:
            #opening the coordinate file
            coordinate_file = open(file, "r")
  
            #reading the file
            coordinate_data = coordinate_file.read()
  
            #coordinates to list
            coordinates = coordinate_data.split("\n")

            coordinate_file.close()

            break
            
        else:
            continue


    #adding column
    df['r6_site'] = pd.Series(coordinates)

    col_list = list(df.columns)

    col_list.insert(1, col_list.pop(-1))

    df = df[col_list]

    print("inserted r6 sites in ", infile)

    #saving as csv
    df.to_csv(infile, index=False)

    print("saved ", infile)


#for parallele mapping
txt_list = ['Chr3L_Zim_RAL_ZI_Fst.csv', 'Chr2L_Zim_RAL_ZI_Fst.csv', 'Chr2R_ZW_RAL_ZI_Fst.csv', 'Chr2L_ZW_RAL_ZI_Fst.csv', 'Chr3L_ZS_RAL_ZI_FR_SAfr_Fst.csv', 'Chr2R_ZS_RAL_ZI_FR_SAfr_Fst.csv', 'ChrX_ZS_ZH_ZW_Fst.csv', 'ChrX_ZS_RAL_ZI_Fst.csv', 'Chr3L_ZH_RAL_ZI_Fst.csv', 'Chr2L_ZS_RAL_ZI_Fst.csv', 'Chr2R_ZS_RAL_ZI_Fst.csv', 'Chr3R_ZH_RAL_ZI_Fst.csv', 'Chr2L_ZS_RAL_ZI_FR_SAfr_Fst.csv', 'Chr3R_ZS_RAL_ZI_FR_SAfr_Fst.csv', 'ChrX_ZW_RAL_ZI_Fst.csv', 'Chr2L_ZS_ZH_ZW_Fst.csv', 'ChrX_Zim_RAL_ZI_Fst.csv', 'Chr3R_ZS_RAL_ZI_Fst.csv', 'Chr3R_Zim_RAL_ZI_Fst.csv', 'Chr2R_ZH_RAL_ZI_Fst.csv', 'Chr2L_ZH_RAL_ZI_Fst.csv', 'Chr3L_ZS_RAL_ZI_Fst.csv', 'Chr2R_Zim_RAL_ZI_Fst.csv', 'Chr3L_ZS_ZH_ZW_Fst.csv', 'Chr2R_ZS_ZH_ZW_Fst.csv', 'ChrX_ZH_RAL_ZI_Fst.csv', 'Chr3R_ZS_ZH_ZW_Fst.csv', 'ChrX_ZS_RAL_ZI_FR_SAfr_Fst.csv', 'Chr3L_ZW_RAL_ZI_Fst.csv', 'Chr3R_ZW_RAL_ZI_Fst.csv']

#run in parallel
def run_in_parallel():
    pool = Pool(processes=8)
    pool.map(coordinate_switcher, txt_list)


if __name__ == '__main__':
    run_in_parallel()


