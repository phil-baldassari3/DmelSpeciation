#script to compute sliding window Fst from per site Fst

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

step = 10000

#sliding window function
def sliding_window_fst(infile):

    if 'X' in infile:
        chrom = 'X'
    elif '2L' in infile:
        chrom = '2L'
    elif '2R' in infile:
        chrom = '2R'
    elif '3L' in infile:
        chrom = '3L'
    elif '3R' in infile:
        chrom = '3R'
    

    #opening dataframe
    df = pd.read_csv(infile)

    print("opened ", infile)

    #sorry, due to the format of my files, this is the only way I can think to implement this script.
    if 'ZS_RAL_ZI_FR_SAfr' not in infile:

        #setting avg fst lists
        avg_fst_1 = []
        avg_fst_2 = []

        #window
        win_start = []
        win_end = []

        print("starting sliding window for ", infile, " with ", step, " bp steps...")

        #sliding window
        ls = df['r6_site'].tolist()
        len_list = [value for value in ls if np.isnan(value) == False]
        print(len_list[-10:-1])


        for i in range(0, int(len_list[-1]+step), step):

            #window lists
            win_start.append(i+1)
            win_end.append(i+step)

            #setting fst wondow list to be averaged
            fst_ls_1 = []
            fst_ls_2 = []

            #looping through dataframe
            for row in df.itertuples():

                row_list = list(row)

                if (row_list[2] >= i+1) & (row_list[2] < i+step):
                    fst_ls_1.append(row_list[3])
                    fst_ls_2.append(row_list[4])
                elif row_list[2] < i+1:
                    continue
                elif row_list[2] > i+step:
                    break
                else:
                    continue

                

            if len(fst_ls_1) == 0:
                avg_fst_1.append(0)
                avg_fst_2.append(0)

            else:
                avg_fst_1.append(sum(fst_ls_1)/len(fst_ls_1))
                avg_fst_2.append(sum(fst_ls_2)/len(fst_ls_2))

        print("making new dataframe for ", infile)

        col_list = list(df.columns)
        win_list = ['window_start', 'window_end']

        
        key_list = win_list + col_list[2:4]
        value_list = [win_start, win_end, avg_fst_1, avg_fst_2]

        dictionary = dict(zip(key_list, value_list))

        windowed_df = pd.DataFrame.from_dict(dictionary)


    else:

        #setting avg fst lists
        avg_fst_1 = []
        avg_fst_2 = []
        avg_fst_3 = []
        avg_fst_4 = []

        #window
        win_start = []
        win_end = []

        print("starting sliding window for ", infile, " with ", step, " bp steps...")

        #sliding window
        ls = df['r6_site'].tolist()
        len_list = [value for value in ls if np.isnan(value) == False]

        for i in range(0, len_list[-1]+step, step):

            #window lists
            win_start.append(i+1)
            win_end.append(i+step)

            #setting fst wondow list to be averaged
            fst_ls_1 = []
            fst_ls_2 = []
            fst_ls_3 = []
            fst_ls_4 = []

            #looping through dataframe
            for row in df.itertuples():

                row_list = list(row)

                if (row_list[2] >= i+1) & (row_list[2] < i+step):
                    fst_ls_1.append(row_list[3])
                    fst_ls_2.append(row_list[4])
                    fst_ls_3.append(row_list[5])
                    fst_ls_4.append(row_list[6])

                elif row_list[2] < i+1:
                    continue
                elif row_list[2] > i+step:
                    break
                else:
                    continue


            if len(fst_ls_1) == 0:
                avg_fst_1.append(0)
                avg_fst_2.append(0)
                avg_fst_3.append(0)
                avg_fst_4.append(0)

            else:
                avg_fst_1.append(sum(fst_ls_1)/len(fst_ls_1))
                avg_fst_2.append(sum(fst_ls_2)/len(fst_ls_2))
                avg_fst_3.append(sum(fst_ls_3)/len(fst_ls_3))
                avg_fst_4.append(sum(fst_ls_4)/len(fst_ls_4))

        print("making new dataframe for ", infile)

        col_list = list(df.columns)
        win_list = ['window_start', 'window_end']

        
        key_list = win_list + col_list[2:6]
        value_list = [win_start, win_end, avg_fst_1, avg_fst_2, avg_fst_3, avg_fst_4]

        dictionary = dict(zip(key_list, value_list))

        windowed_df = pd.DataFrame.from_dict(dictionary)

    #adding chrom column
    windowed_df.insert(0, 'Chrom', chrom)

    windowed_df.to_csv('windowed_10kbp_{file}'.format(file=infile), index=False)

    print("saved dataframe.")


    print("DO AVERAGES IN R. I DON'T FEEL LIKE SCRIPTING IT HERE.")



#for parallel mapping
csv_list = ['Chr2R_ZW_RAL_ZI_Fst.csv', 'ChrX_ZS_ZH_ZW_Fst.csv', 'Chr2R_ZH_RAL_ZI_Fst.csv', 'Chr2R_Zim_RAL_ZI_Fst.csv']

#run in parallel
def run_in_parallel():
    pool = Pool(processes=4)
    pool.map(sliding_window_fst, csv_list)


if __name__ == '__main__':
    run_in_parallel()




    
