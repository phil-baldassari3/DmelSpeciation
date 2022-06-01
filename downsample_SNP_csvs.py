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

    #set Chrom
    chrom = 'Chr'
    
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
    loci.drop(loci.columns.difference(['Locus']), axis=1, inplace=True)

    #open dataframe
    df = pd.read_csv(infile)

    print('looping through rows and downsampling by population for ' + infile)

    #loop through rows
    for row in df.itertuples():

        #tuple to list and counts
        row_list = list(row)

        Acount = row_list.count("A") 
        Tcount = row_list.count("T")
        Ccount = row_list.count("C")
        Gcount = row_list.count("G")

        if Acount >= Tcount and Acount >= Ccount and Acount >= Gcount:
            major_allele = "A"
        elif Tcount >= Acount and Tcount >= Ccount and Tcount >= Gcount:
            major_allele = "T"
        elif Ccount >= Acount and Ccount >= Tcount and Ccount >= Gcount:
            major_allele = "C"
        elif Gcount >= Acount and Gcount >= Tcount and Gcount >= Ccount:
            major_allele = "G"

        ###list for each population

        ##FR

        #list downsampling
        FR_row = row_list[1:88]
        no_N_FR_row = [value for value in FR_row if value != 'N']
        random_FR_row = random.sample(no_N_FR_row, len(no_N_FR_row))
        downsample_FR = random_FR_row[0:FR_n]

        #sample size
        bpcount_FR = len(downsample_FR)

        #maf
        if bpcount_FR == 0:
            freq_FR = np.nan
        else:
            freq_FR = 1 - (downsample_FR.count(major_allele) / bpcount_FR)

        #appending
        n_FR.append(bpcount_FR)
        maf_FR.append(freq_FR)

        ##RAL

        #list downsampling
        RAL_row = row_list[88:293]
        no_N_RAL_row = [value for value in RAL_row if value != 'N']
        random_RAL_row = random.sample(no_N_RAL_row, len(no_N_RAL_row))
        downsample_RAL = random_RAL_row[0:RAL_n]

        #sample size
        bpcount_RAL = len(downsample_RAL)

        #maf
        if bpcount_RAL == 0:
            freq_RAL = np.nan
        else:
            freq_RAL = 1 - (downsample_RAL.count(major_allele) / bpcount_RAL)

        #appending
        n_RAL.append(bpcount_RAL)
        maf_RAL.append(freq_RAL)

        ##SAfr

        #list downsampling
        SAfr_row = row_list[293:406]
        no_N_SAfr_row = [value for value in SAfr_row if value != 'N']
        random_SAfr_row = random.sample(no_N_SAfr_row, len(no_N_SAfr_row))
        downsample_SAfr = random_SAfr_row[0:SAfr_n]

        #sample size
        bpcount_SAfr = len(downsample_SAfr)

        #maf
        if bpcount_SAfr == 0:
            freq_SAfr = np.nan
        else:
            freq_SAfr = 1 - (downsample_SAfr.count(major_allele) / bpcount_SAfr)

        #appending
        n_SAfr.append(bpcount_SAfr)
        maf_SAfr.append(freq_SAfr)

        ##ZI

        #list downsampling
        ZI_row = row_list[410:606]
        no_N_ZI_row = [value for value in ZI_row if value != 'N']
        random_ZI_row = random.sample(no_N_ZI_row, len(no_N_ZI_row))
        downsample_ZI = random_ZI_row[0:ZI_n]

        #sample size
        bpcount_ZI = len(downsample_ZI)

        #maf
        if bpcount_ZI == 0:
            freq_ZI = np.nan
        else:
            freq_ZI = 1 - (downsample_ZI.count(major_allele) / bpcount_ZI)

        #appending
        n_ZI.append(bpcount_ZI)
        maf_ZI.append(freq_ZI)

        ##ZH

        #list downsampling
        ZH_row = row_list[406:410]
        no_N_ZH_row = [value for value in ZH_row if value != 'N']
        random_ZH_row = random.sample(no_N_ZH_row, len(no_N_ZH_row))
        downsample_ZH = random_ZH_row[0:ZH_n]

        #sample size
        bpcount_ZH = len(downsample_ZH)

        #maf
        if bpcount_ZH == 0:
            freq_ZH = np.nan
        else:
            freq_ZH = 1 - (downsample_ZH.count(major_allele) / bpcount_ZH)

        #appending
        n_ZH.append(bpcount_ZH)
        maf_ZH.append(freq_ZH)

        ##ZW

        #list downsampling
        ZW_row = row_list[611:620]
        no_N_ZW_row = [value for value in ZW_row if value != 'N']
        random_ZW_row = random.sample(no_N_ZW_row, len(no_N_ZW_row))
        downsample_ZW = random_ZW_row[0:ZW_n]

        #sample size
        bpcount_ZW = len(downsample_ZW)

        #maf
        if bpcount_ZW == 0:
            freq_ZW = np.nan
        else:
            freq_ZW = 1 - (downsample_ZW.count(major_allele) / bpcount_ZW)

        #appending
        n_ZW.append(bpcount_ZW)
        maf_ZW.append(freq_ZW)

        ##ZS

        #list downsampling
        ZS_row = row_list[606:611]
        no_N_ZS_row = [value for value in ZS_row if value != 'N']
        random_ZS_row = random.sample(no_N_ZS_row, len(no_N_ZS_row))
        downsample_ZS = random_ZS_row[0:ZS_n]

        #sample size
        bpcount_ZS = len(downsample_ZS)

        #maf
        if bpcount_ZS == 0:
            freq_ZS = np.nan
        else:
            freq_ZS = 1 - (downsample_ZS.count(major_allele) / bpcount_ZS)

        #appending
        n_ZS.append(bpcount_ZS)
        maf_ZS.append(freq_ZS)

    #dictionary to dataframe
    dictionary = {'n_FR': n_FR, 'maf_FR': maf_FR, 'n_RAL': n_RAL, 'maf_RAL': maf_RAL, 'n_SAfr': n_SAfr, 'maf_SAfr': maf_SAfr, 'n_ZI': n_ZI, 'maf_ZI': maf_ZI, 'n_ZH': n_ZH, 'maf_ZH': maf_ZH, 'n_ZW': n_ZW, 'maf_ZW': maf_ZW, 'n_ZS': n_ZS, 'maf_ZS': maf_ZS}

    data_df = pd.DataFrame(dictionary)

    #merge dataframes
    n_freq_df = pd.merge(loci, data_df, left_index=True, right_index=True)

    #adding chrom column
    n_freq_df.insert(0, 'Chrom', chrom)

    #zim composite populations
    n_freq_df['n_ZH_ZW'] = n_freq_df['n_ZH'] + n_freq_df['n_ZW']
    n_freq_df['maf_ZH_ZW'] = ((n_freq_df['maf_ZH'] * n_freq_df['n_ZH']) + (n_freq_df['maf_ZW'] * n_freq_df['n_ZW'])) / n_freq_df['n_ZH_ZW']
    n_freq_df['n_zim'] = n_freq_df['n_ZH'] + n_freq_df['n_ZW'] + n_freq_df['n_ZS']
    n_freq_df['maf_zim'] = ((n_freq_df['maf_ZH'] * n_freq_df['n_ZH']) + (n_freq_df['maf_ZW'] * n_freq_df['n_ZW']) * (n_freq_df['maf_ZS'] * n_freq_df['n_ZS'])) / n_freq_df['n_zim']


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





