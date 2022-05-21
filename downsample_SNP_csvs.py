#script to downsample SNP data from the SNP csv files. will return dataframes with allele frequencies and population counts. data must be filtered first.

#importing modules
import os,sys
import pandas as pd
import numpy as np
from multiprocessing import Pool
import random

#setting directory
directory = "set"

'''
#get list of csvs for parallel pool mapping
file_list = []

for file in os.listdir(directory):
    if file.endswith('.csv'):
        file_list.append(file)
    else:
        continue

print(file_list)
'''

#defining downsampled n for each pop
FR_n = 40
RAL_n = 100
SAfr_n = 40
ZI_n = 100
ZH_n = 3
ZW_n = 6
ZS_n = 4

#downsampling function
def snp_csv_downsample(infile):

    print('starting process for ' + infile)
    
    #lists
    n_FR = []
    n_RAL = []
    n_SAfr = []
    n_ZI = []
    n_ZH = []
    n_ZW = []
    n_ZS = []
    maf_FR = []
    maf_RAL = []
    maf_SAfr = []
    maf_ZI = []
    maf_ZH = []
    maf_ZW = []
    maf_ZS = []

    #open dataframe
    loci = pd.read_csv(infile)
    #locus dataframe
    loci = loci.drop(loci.columns.difference(['Locus']), axis=1, inplace=True)

    #open dataframe
    df = pd.read_csv(infile)
    #sequence dataframe
    df = df.drop(['Locus', 'Acount', 'Tcount', 'Ccount', 'Gcount', 'bpcount', 'Aprop', 'Tprop', 'Cprop', 'Gprop', 'MajorAF'], axis=1)

    print('looping through rows and downsampling by population for ' + infile)

    #loop through rows
    for row in df.itertuples():

        #tuple to list
        row_list = list(row)

        ###list for each population

        ##FR

        #list downsampling
        FR_row = row_list[0:87]
        FR_row.append('X')
        no_N_FR_row = FR_row.remove('N')
        random_FR_row = random.sample(no_N_FR_row, len(no_N_FR_row))
        random_FR_row.append('Y')
        random_FR_row.remove('X')
        downsample_FR = random_FR_row[0:FR_n]

        #sample size
        bpcount_FR = len(downsample_FR)

        #count alleles, maf
        Acount_FR = downsample_FR.count("A")
        Tcount_FR = downsample_FR.count("T")
        Ccount_FR = downsample_FR.count("C")
        Gcount_FR = downsample_FR.count("G")

        count_list = [Acount_FR, Tcount_FR, Ccount_FR, Gcount_FR]

        freq_FR = 1 - (max(count_list) / bpcount_FR)

        #appending
        n_FR.append(bpcount_FR)
        maf_FR.append(freq_FR)

        ##RAL

        #list downsampling
        RAL_row = row_list[87:292]
        RAL_row.append('X')
        no_N_RAL_row = RAL_row.remove('N')
        random_RAL_row = random.sample(no_N_RAL_row, len(no_N_RAL_row))
        random_RAL_row.append('Y')
        random_RAL_row.remove('X')
        downsample_RAL = random_RAL_row[0:RAL_n]

        #sample size
        bpcount_RAL = len(downsample_RAL)

        #count alleles, maf
        Acount_RAL = downsample_RAL.count("A")
        Tcount_RAL = downsample_RAL.count("T")
        Ccount_RAL = downsample_RAL.count("C")
        Gcount_RAL = downsample_RAL.count("G")

        count_list = [Acount_RAL, Tcount_RAL, Ccount_RAL, Gcount_RAL]

        freq_RAL = 1 - (max(count_list) / bpcount_RAL)

        #appending
        n_RAL.append(bpcount_RAL)
        maf_RAL.append(freq_RAL)

        ##SAfr

        #list downsampling
        SAfr_row = row_list[292:405]
        SAfr_row.append('X')
        no_N_SAfr_row = SAfr_row.remove('N')
        random_SAfr_row = random.sample(no_N_SAfr_row, len(no_N_SAfr_row))
        random_SAfr_row.append('Y')
        random_SAfr_row.remove('X')
        downsample_SAfr = random_SAfr_row[0:SAfr_n]

        #sample size
        bpcount_SAfr = len(downsample_SAfr)

        #count alleles, maf
        Acount_SAfr = downsample_SAfr.count("A")
        Tcount_SAfr = downsample_SAfr.count("T")
        Ccount_SAfr = downsample_SAfr.count("C")
        Gcount_SAfr = downsample_SAfr.count("G")

        count_list = [Acount_SAfr, Tcount_SAfr, Ccount_SAfr, Gcount_SAfr]

        freq_SAfr = 1 - (max(count_list) / bpcount_SAfr)

        #appending
        n_SAfr.append(bpcount_SAfr)
        maf_SAfr.append(freq_SAfr)

        ##ZI

        #list downsampling
        ZI_row = row_list[409:605]
        ZI_row.append('X')
        no_N_ZI_row = ZI_row.remove('N')
        random_ZI_row = random.sample(no_N_ZI_row, len(no_N_ZI_row))
        random_ZI_row.append('Y')
        random_ZI_row.remove('X')
        downsample_ZI = random_ZI_row[0:ZI_n]

        #sample size
        bpcount_ZI = len(downsample_ZI)

        #count alleles, maf
        Acount_ZI = downsample_ZI.count("A")
        Tcount_ZI = downsample_ZI.count("T")
        Ccount_ZI = downsample_ZI.count("C")
        Gcount_ZI = downsample_ZI.count("G")

        count_list = [Acount_ZI, Tcount_ZI, Ccount_ZI, Gcount_ZI]

        freq_ZI = 1 - (max(count_list) / bpcount_ZI)

        #appending
        n_ZI.append(bpcount_ZI)
        maf_ZI.append(freq_ZI)

        ##ZH

        #list downsampling
        ZH_row = row_list[405:409]
        ZH_row.append('X')
        no_N_ZH_row = ZH_row.remove('N')
        random_ZH_row = random.sample(no_N_ZH_row, len(no_N_ZH_row))
        random_ZH_row.append('Y')
        random_ZH_row.remove('X')
        downsample_ZH = random_ZH_row[0:ZH_n]

        #sample size
        bpcount_ZH = len(downsample_ZH)

        #count alleles, maf
        Acount_ZH = downsample_ZH.count("A")
        Tcount_ZH = downsample_ZH.count("T")
        Ccount_ZH = downsample_ZH.count("C")
        Gcount_ZH = downsample_ZH.count("G")

        count_list = [Acount_ZH, Tcount_ZH, Ccount_ZH, Gcount_ZH]

        freq_ZH = 1 - (max(count_list) / bpcount_ZH)

        #appending
        n_ZH.append(bpcount_ZH)
        maf_ZH.append(freq_ZH)

        ##ZW

        #list downsampling
        ZW_row = row_list[610:619]
        ZW_row.append('X')
        no_N_ZW_row = ZW_row.remove('N')
        random_ZW_row = random.sample(no_N_ZW_row, len(no_N_ZW_row))
        random_ZW_row.append('Y')
        random_ZW_row.remove('X')
        downsample_ZW = random_ZW_row[0:ZW_n]

        #sample size
        bpcount_ZW = len(downsample_ZW)

        #count alleles, maf
        Acount_ZW = downsample_ZW.count("A")
        Tcount_ZW = downsample_ZW.count("T")
        Ccount_ZW = downsample_ZW.count("C")
        Gcount_ZW = downsample_ZW.count("G")

        count_list = [Acount_ZW, Tcount_ZW, Ccount_ZW, Gcount_ZW]

        freq_ZW = 1 - (max(count_list) / bpcount_ZW)

        #appending
        n_ZW.append(bpcount_ZW)
        maf_ZW.append(freq_ZW)

        ##ZS

        #list downsampling
        ZS_row = row_list[605:610]
        ZS_row.append('X')
        no_N_ZS_row = ZS_row.remove('N')
        random_ZS_row = random.sample(no_N_ZS_row, len(no_N_ZS_row))
        random_ZS_row.append('Y')
        random_ZS_row.remove('X')
        downsample_ZS = random_ZS_row[0:ZS_n]

        #sample size
        bpcount_ZS = len(downsample_ZS)

        #count alleles, maf
        Acount_ZS = downsample_ZS.count("A")
        Tcount_ZS = downsample_ZS.count("T")
        Ccount_ZS = downsample_ZS.count("C")
        Gcount_ZS = downsample_ZS.count("G")

        count_list = [Acount_ZS, Tcount_ZS, Ccount_ZS, Gcount_ZS]

        freq_ZS = 1 - (max(count_list) / bpcount_ZS)

        #appending
        n_ZS.append(bpcount_ZS)
        maf_ZS.append(freq_ZS)

    #dictionary to dataframe
    dictionary = {'n_FR': n_FR, 'maf_FR': maf_FR, 'n_RAL': n_RAL, 'maf_RAL': maf_RAL, 'n_SAfr': n_SAfr, 'maf_SAfr': maf_SAfr, 'n_ZI': n_ZI, 'maf_ZI': maf_ZI, 'n_ZH': n_ZH, 'maf_ZH': maf_ZH, 'n_ZW': n_ZW, 'maf_ZW': maf_ZW, 'n_ZS': n_ZS, 'maf_ZS': maf_ZS}

    data_df = pd.DataFrame(dictionary)

    #merge dataframes
    n_freq_df = pd.merge(loci, data_df, left_index=True, right_index=True)

    print('new dataframe made for ' + infile)

    #create csv file
    number = infile.split("_")[3]

    n_freq_df.to_csv('downsampled_allele_freq_chunk{num}.csv'.format(num=number), index=False)

    print('new file made from ' + infile)




#for parallele mapping
csv_list = []

#run in parallel
def run_in_parallel():
    pool = Pool(processes=8)
    pool.map(snp_csv_downsample, csv_list)


if __name__ == '__main__':
    run_in_parallel()





