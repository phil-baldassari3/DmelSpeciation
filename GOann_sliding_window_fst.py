#script to compute sliding window Fst from per site Fst

#importing modules
import os,sys
import pandas as pd
import numpy as np
from multiprocessing import Pool

#setting directory
directory = "/Users/philipbaldassari/Desktop/zim-cos_gene_GO_ann"
os.chdir(directory)

'''
#get list of txt list files for parallel pool mapping
file_list = []

for file in os.listdir(directory):
    if file.endswith('.csv') and 'GO_sites_genes' in file:
        file_list.append(file)
    else:
        continue

print(file_list)
'''

step = 1000



def sliding_window(infile, df):

    """Takes annotated df of per site Fst and returns a sliding window Fst annotated. If this is an autosome dataframe you need to run it through autosome_splitter() first"""

    #setting chrom
    chrom_ls = df['Chrom'].tolist()

    if 'ChrX' in chrom_ls:
        chrom = 'X'
    elif 'Chr2L' in chrom_ls:
        chrom = '2L'
    elif 'Chr2R' in chrom_ls:
        chrom = '2R'
    elif 'Chr3L' in chrom_ls:
        chrom = '3L'
    elif 'Chr3R' in chrom_ls:
        chrom = '3R'
    

    

    #sorry, due to the format of my files, this is the only way I can think to implement this script.
    if 'ZS_RAL_ZI_FR_SAfr' in infile:

        #setting avg fst lists
        avg_fst_1 = []
        avg_fst_2 = []
        avg_fst_3 = []
        avg_fst_4 = []

        #window
        win_start = []
        win_end = []

        #setting annotation lists
        FBgn = []
        gene_symbol = []
        neurogenesis = []
        mating_behavior = []
        male_mating_behavior = []
        female_mating_behavior = []


        print("starting sliding window for ", infile, " with ", step, " bp steps...")

        #sliding window
        len_list = df['r6_site'].tolist()

        for i in range(0, int(len_list[-1])+step, step):

           #window lists
            win_start.append(i+1)
            win_end.append(i+step)

            #chunking df
            chunk_df = df.drop(df[~((df.r6_site >= i+1) & (df.r6_site < i+step))].index)

            if len(chunk_df.index) == 0:
                avg_fst_1.append(0)
                avg_fst_2.append(0)
                avg_fst_3.append(0)
                avg_fst_4.append(0)

                FBgn.append('None')
                gene_symbol.append('None')
                neurogenesis.append('None')
                mating_behavior.append('None')
                male_mating_behavior.append('None')
                female_mating_behavior.append('None')

            else:
                avg_fst_1.append(chunk_df.iloc[:, 2].mean())
                avg_fst_2.append(chunk_df.iloc[:, 3].mean())
                avg_fst_3.append(chunk_df.iloc[:, 4].mean())
                avg_fst_4.append(chunk_df.iloc[:, 5].mean())

                chunk_FBgn_temp = chunk_df['FBgn'].tolist()
                str_FBgn = "; ".join(chunk_FBgn_temp)
                chunk_FBgn = list(set(str_FBgn.split('; ')))
                if len(chunk_FBgn) == 1:
                    FBgn.append(''.join(chunk_FBgn))
                else:
                    chunk_FBgn = [i for i in chunk_FBgn if i != "None"]
                    FBgn.append("; ".join(chunk_FBgn))
                
                chunk_gene_symbol_temp = chunk_df['gene_symbol'].tolist()
                chunk_gene_symbol_temp = [str(i) for i in chunk_gene_symbol_temp]
                str_gene_symbol = "; ".join(chunk_gene_symbol_temp)
                chunk_gene_symbol = list(set(str_gene_symbol.split('; ')))
                if len(chunk_gene_symbol) == 1:
                    gene_symbol.append(''.join(chunk_gene_symbol))
                else:
                    chunk_gene_symbol = [i for i in chunk_gene_symbol if i != "None"]
                    gene_symbol.append("; ".join(chunk_gene_symbol))
                
                chunk_neurogenesis_temp = chunk_df['Neurogenesis'].tolist()
                str_neurogenesis = "; ".join(chunk_neurogenesis_temp)
                chunk_neurogenesis = list(set(str_neurogenesis.split('; ')))
                if 'Neurogenesis' in chunk_neurogenesis:
                    neurogenesis.append('Neurogenesis')
                else:
                    neurogenesis.append('None')


                chunk_mating_behavior_temp = chunk_df['Mating_behavior'].tolist()
                str_mating_behavior = "; ".join(chunk_mating_behavior_temp)
                chunk_mating_behavior = list(set(str_mating_behavior.split('; ')))
                if 'Mating_behavior' in chunk_mating_behavior:
                    mating_behavior.append('Mating_behavior')
                else:
                    mating_behavior.append('None')

                chunk_male_mating_behavior_temp = chunk_df['Male_mating_behavior'].tolist()
                str_male_mating_behavior = "; ".join(chunk_male_mating_behavior_temp)
                chunk_male_mating_behavior = list(set(str_male_mating_behavior.split('; ')))
                if 'Male_mating_behavior' in chunk_male_mating_behavior:
                    male_mating_behavior.append('Male_mating_behavior')
                else:
                    male_mating_behavior.append('None')

                chunk_female_mating_behavior_temp = chunk_df['Female_mating_behavior'].tolist()
                str_female_mating_behavior = "; ".join(chunk_female_mating_behavior_temp)
                chunk_female_mating_behavior = list(set(str_female_mating_behavior.split('; ')))
                if 'Female_mating_behavior' in chunk_female_mating_behavior:
                    female_mating_behavior.append('Female_mating_behavior')
                else:
                    female_mating_behavior.append('None')


        print("making new dataframe for ", infile)

        #key lists
        col_list = list(df.columns)
        win_list = ['window_start', 'window_end']
        ann_list = ['FBgn', 'gene_symbol', 'Neurogenesis', 'Mating_behavior', 'Male_mating_behavior', 'Female_mating_behavior']

        
        key_list = win_list + col_list[2:6] + ann_list

        #value list
        value_list = [win_start, win_end, avg_fst_1, avg_fst_2, avg_fst_3, avg_fst_4, FBgn, gene_symbol, neurogenesis, mating_behavior, male_mating_behavior, female_mating_behavior]

        #merging key, values to dictionary
        dictionary = dict(zip(key_list, value_list))

        #making df
        windowed_df = pd.DataFrame.from_dict(dictionary)

  


    elif 'RAL_vs_FR_ZI_SAfr' in infile or 'RAL_vs_FR_ZI_SAfr' in infile or 'RAL_vs_FR_ZI_SAfr' in infile or 'RAL_vs_FR_ZI_SAfr' in infile:

        #setting avg fst lists
        avg_fst_1 = []
        avg_fst_2 = []
        avg_fst_3 = []

        #window
        win_start = []
        win_end = []

        #setting annotation lists
        FBgn = []
        gene_symbol = []
        neurogenesis = []
        mating_behavior = []
        male_mating_behavior = []
        female_mating_behavior = []

        print("starting sliding window for ", infile, " with ", step, " bp steps...")

        #sliding window
        len_list = df['r6_site'].tolist()

        for i in range(0, int(len_list[-1])+step, step):

           #window lists
            win_start.append(i+1)
            win_end.append(i+step)

            #chunking df
            chunk_df = df.drop(df[~((df.r6_site >= i+1) & (df.r6_site < i+step))].index)

            if len(chunk_df.index) == 0:
                avg_fst_1.append(0)
                avg_fst_2.append(0)
                avg_fst_3.append(0)

                FBgn.append('None')
                gene_symbol.append('None')
                neurogenesis.append('None')
                mating_behavior.append('None')
                male_mating_behavior.append('None')
                female_mating_behavior.append('None')

            else:
                avg_fst_1.append(chunk_df.iloc[:, 2].mean())
                avg_fst_2.append(chunk_df.iloc[:, 3].mean())
                avg_fst_3.append(chunk_df.iloc[:, 4].mean())

                chunk_FBgn_temp = chunk_df['FBgn'].tolist()
                str_FBgn = "; ".join(chunk_FBgn_temp)
                chunk_FBgn = list(set(str_FBgn.split('; ')))
                if len(chunk_FBgn) == 1:
                    FBgn.append(''.join(chunk_FBgn))
                else:
                    chunk_FBgn = [i for i in chunk_FBgn if i != "None"]
                    FBgn.append("; ".join(chunk_FBgn))
                
                chunk_gene_symbol_temp = chunk_df['gene_symbol'].tolist()
                chunk_gene_symbol_temp = [str(i) for i in chunk_gene_symbol_temp]
                str_gene_symbol = "; ".join(chunk_gene_symbol_temp)
                chunk_gene_symbol = list(set(str_gene_symbol.split('; ')))
                if len(chunk_gene_symbol) == 1:
                    gene_symbol.append(''.join(chunk_gene_symbol))
                else:
                    chunk_gene_symbol = [i for i in chunk_gene_symbol if i != "None"]
                    gene_symbol.append("; ".join(chunk_gene_symbol))
                
                chunk_neurogenesis_temp = chunk_df['Neurogenesis'].tolist()
                str_neurogenesis = "; ".join(chunk_neurogenesis_temp)
                chunk_neurogenesis = list(set(str_neurogenesis.split('; ')))
                if 'Neurogenesis' in chunk_neurogenesis:
                    neurogenesis.append('Neurogenesis')
                else:
                    neurogenesis.append('None')


                chunk_mating_behavior_temp = chunk_df['Mating_behavior'].tolist()
                str_mating_behavior = "; ".join(chunk_mating_behavior_temp)
                chunk_mating_behavior = list(set(str_mating_behavior.split('; ')))
                if 'Mating_behavior' in chunk_mating_behavior:
                    mating_behavior.append('Mating_behavior')
                else:
                    mating_behavior.append('None')

                chunk_male_mating_behavior_temp = chunk_df['Male_mating_behavior'].tolist()
                str_male_mating_behavior = "; ".join(chunk_male_mating_behavior_temp)
                chunk_male_mating_behavior = list(set(str_male_mating_behavior.split('; ')))
                if 'Male_mating_behavior' in chunk_male_mating_behavior:
                    male_mating_behavior.append('Male_mating_behavior')
                else:
                    male_mating_behavior.append('None')

                chunk_female_mating_behavior_temp = chunk_df['Female_mating_behavior'].tolist()
                str_female_mating_behavior = "; ".join(chunk_female_mating_behavior_temp)
                chunk_female_mating_behavior = list(set(str_female_mating_behavior.split('; ')))
                if 'Female_mating_behavior' in chunk_female_mating_behavior:
                    female_mating_behavior.append('Female_mating_behavior')
                else:
                    female_mating_behavior.append('None')


        print("making new dataframe for ", infile)

        col_list = list(df.columns)
        win_list = ['window_start', 'window_end']
        ann_list = ['FBgn', 'gene_symbol', 'Neurogenesis', 'Mating_behavior', 'Male_mating_behavior', 'Female_mating_behavior']
        
        key_list = win_list + col_list[2:5] + ann_list
        value_list = [win_start, win_end, avg_fst_1, avg_fst_2, avg_fst_3, FBgn, gene_symbol, neurogenesis, mating_behavior, male_mating_behavior, female_mating_behavior]
        

        dictionary = dict(zip(key_list, value_list))

        windowed_df = pd.DataFrame.from_dict(dictionary)



    else:

        #setting avg fst lists
        avg_fst_1 = []
        avg_fst_2 = []

        #window
        win_start = []
        win_end = []

        #setting annotation lists
        FBgn = []
        gene_symbol = []
        neurogenesis = []
        mating_behavior = []
        male_mating_behavior = []
        female_mating_behavior = []

        print("starting sliding window for ", infile, " with ", step, " bp steps...")

        #sliding window
        len_list = df['r6_site'].tolist()
        
        for i in range(0, int(len_list[-1])+step, step):

            #window lists
            win_start.append(i+1)
            win_end.append(i+step)

            #chunking df
            chunk_df = df.drop(df[~((df.r6_site >= i+1) & (df.r6_site < i+step))].index)

            if len(chunk_df.index) == 0:
                avg_fst_1.append(0)
                avg_fst_2.append(0)

                FBgn.append('None')
                gene_symbol.append('None')
                neurogenesis.append('None')
                mating_behavior.append('None')
                male_mating_behavior.append('None')
                female_mating_behavior.append('None')

            else:
                avg_fst_1.append(chunk_df.iloc[:, 2].mean())
                avg_fst_2.append(chunk_df.iloc[:, 3].mean())

                chunk_FBgn_temp = chunk_df['FBgn'].tolist()
                str_FBgn = "; ".join(chunk_FBgn_temp)
                chunk_FBgn = list(set(str_FBgn.split('; ')))
                if len(chunk_FBgn) == 1:
                    FBgn.append(''.join(chunk_FBgn))
                else:
                    chunk_FBgn = [i for i in chunk_FBgn if i != "None"]
                    FBgn.append("; ".join(chunk_FBgn))
                
                chunk_gene_symbol_temp = chunk_df['gene_symbol'].tolist()
                chunk_gene_symbol_temp = [str(i) for i in chunk_gene_symbol_temp]
                str_gene_symbol = "; ".join(chunk_gene_symbol_temp)
                chunk_gene_symbol = list(set(str_gene_symbol.split('; ')))
                if len(chunk_gene_symbol) == 1:
                    gene_symbol.append(''.join(chunk_gene_symbol))
                else:
                    chunk_gene_symbol = [i for i in chunk_gene_symbol if i != "None"]
                    gene_symbol.append("; ".join(chunk_gene_symbol))
                
                chunk_neurogenesis_temp = chunk_df['Neurogenesis'].tolist()
                str_neurogenesis = "; ".join(chunk_neurogenesis_temp)
                chunk_neurogenesis = list(set(str_neurogenesis.split('; ')))
                if 'Neurogenesis' in chunk_neurogenesis:
                    neurogenesis.append('Neurogenesis')
                else:
                    neurogenesis.append('None')


                chunk_mating_behavior_temp = chunk_df['Mating_behavior'].tolist()
                str_mating_behavior = "; ".join(chunk_mating_behavior_temp)
                chunk_mating_behavior = list(set(str_mating_behavior.split('; ')))
                if 'Mating_behavior' in chunk_mating_behavior:
                    mating_behavior.append('Mating_behavior')
                else:
                    mating_behavior.append('None')

                chunk_male_mating_behavior_temp = chunk_df['Male_mating_behavior'].tolist()
                str_male_mating_behavior = "; ".join(chunk_male_mating_behavior_temp)
                chunk_male_mating_behavior = list(set(str_male_mating_behavior.split('; ')))
                if 'Male_mating_behavior' in chunk_male_mating_behavior:
                    male_mating_behavior.append('Male_mating_behavior')
                else:
                    male_mating_behavior.append('None')

                chunk_female_mating_behavior_temp = chunk_df['Female_mating_behavior'].tolist()
                str_female_mating_behavior = "; ".join(chunk_female_mating_behavior_temp)
                chunk_female_mating_behavior = list(set(str_female_mating_behavior.split('; ')))
                if 'Female_mating_behavior' in chunk_female_mating_behavior:
                    female_mating_behavior.append('Female_mating_behavior')
                else:
                    female_mating_behavior.append('None')

        print("making new dataframe for ", infile)

        col_list = list(df.columns)
        win_list = ['window_start', 'window_end']
        ann_list = ['FBgn', 'gene_symbol', 'Neurogenesis', 'Mating_behavior', 'Male_mating_behavior', 'Female_mating_behavior']

        
        key_list = win_list + col_list[2:4] + ann_list
        value_list = [win_start, win_end, avg_fst_1, avg_fst_2, FBgn, gene_symbol, neurogenesis, mating_behavior, male_mating_behavior, female_mating_behavior]


        """ print('cols:', key_list)
        print('win_start:', win_start[1328:1348], len(win_start))
        print('win_end:', win_end[1328:1348], len(win_end))
        print('avg_fst_1:', avg_fst_1[1328:1348], len(avg_fst_1))
        print('avg_fst_2:', avg_fst_2[1328:1348], len(avg_fst_2))
        print('FBgn:', FBgn[1328:1348], len(FBgn))
        print('gene_symbol:', gene_symbol[1328:1348], len(gene_symbol))
        print('neurogenesis:', neurogenesis[1328:1348], len(neurogenesis))
        print('mating_behavior:', mating_behavior[1328:1348], len(mating_behavior))
        print('male_mating_behavior:', male_mating_behavior[1328:1348], len(male_mating_behavior))
        print('female_mating_behavior:', female_mating_behavior[1328:1348], len(female_mating_behavior)) """



        dictionary = dict(zip(key_list, value_list))

        windowed_df = pd.DataFrame.from_dict(dictionary)




    
    #adding chrom column
    windowed_df.insert(0, 'Chrom', chrom)

    return windowed_df


def autosome_splitter(infile, df):

    """splits an autosome df into searate chromosome dataframes"""
    if "ZH_RAL_ZI" in infile:
        df_2R = df[df.Chrom == "Chr2R"]
        df_3L = df[df.Chrom == "Chr3L"]
        df_3R = df[df.Chrom == "Chr3R"]

        return df_2R, df_3L, df_3R

    else:
        df_2L = df[df.Chrom == "Chr2L"]
        df_2R = df[df.Chrom == "Chr2R"]
        df_3L = df[df.Chrom == "Chr3L"]
        df_3R = df[df.Chrom == "Chr3R"]

        return df_2L, df_2R, df_3L, df_3R


def main_function(infile):

    """main function of the program. uses sliding_window() and autosome_splitter() to output annotated sliding window Fst files for each infile. THis function in run in parallel for each infile"""

    persite_df = pd.read_csv(infile)

    if "X" in infile:
        win_df = sliding_window(infile, persite_df)

        win_df.to_csv('windowed_1kbp_{file}'.format(file=infile), index=False)

        print("saved dataframe.")


        print("DO AVERAGES IN R. I DON'T FEEL LIKE SCRIPTING IT HERE.")



    elif "ZH_RAL_ZI" in infile:
        persite_df_2R, persite_df_3L, persite_df_3R = autosome_splitter(infile, persite_df)

        win_df_2R = sliding_window(infile, persite_df_2R)
        win_df_3L = sliding_window(infile, persite_df_3L)
        win_df_3R = sliding_window(infile, persite_df_3R)

        win_df = pd.concat([win_df_2R, win_df_3L, win_df_3R])

        win_df.to_csv('windowed_1kbp_{file}'.format(file=infile), index=False)

        print("saved dataframe.")
        print("DO AVERAGES IN R. I DON'T FEEL LIKE SCRIPTING IT HERE.")

    else:
        persite_df_2L, persite_df_2R, persite_df_3L, persite_df_3R = autosome_splitter(infile, persite_df)

        win_df_2L = sliding_window(infile, persite_df_2L)
        win_df_2R = sliding_window(infile, persite_df_2R)
        win_df_3L = sliding_window(infile, persite_df_3L)
        win_df_3R = sliding_window(infile, persite_df_3R)

        win_df = pd.concat([win_df_2L, win_df_2R, win_df_3L, win_df_3R])

        win_df.to_csv('windowed_1kbp_{file}'.format(file=infile), index=False)

        print("saved dataframe.")
        print("DO AVERAGES IN R. I DON'T FEEL LIKE SCRIPTING IT HERE.")




#main_function('GO_sites_genes_ZS_RAL_ZI_Fst_autosomes.csv')





#for parallel mapping
csv_list = ['GO_sites_genes_SAfr_vs_RAL_FR_ZI_Fst_autosomes.csv', 'GO_sites_genes_ZS_ZH_ZW_Fst_ChrX.csv', 'GO_sites_genes_FR_vs_RAL_ZI_SAfr_Fst_autosomes.csv', 'GO_sites_genes_FR_vs_RAL_ZI_SAfr_Fst_ChrX.csv', 'GO_sites_genes_ZW_RAL_ZI_Fst_autosomes.csv', 'GO_sites_genes_ZS_RAL_ZI_FR_SAfr_Fst_ChrX.csv', 'GO_sites_genes_ZS_RAL_ZI_Fst_autosomes.csv', 'GO_sites_genes_ZS_RAL_ZI_FR_SAfr_Fst_autosomes.csv', 'GO_sites_genes_ZI_vs_RAL_FR_SAfr_Fst_autosomes.csv', 'GO_sites_genes_ZI_vs_RAL_FR_SAfr_Fst_ChrX.csv', 'GO_sites_genes_ZH_RAL_ZI_Fst_ChrX.csv', 'GO_sites_genes_ZW_RAL_ZI_Fst_ChrX.csv', 'GO_sites_genes_SAfr_vs_RAL_FR_ZI_Fst_ChrX.csv', 'GO_sites_genes_ZH_RAL_ZI_Fst_autosomes.csv', 'GO_sites_genes_ZS_ZH_ZW_Fst_autosomes.csv', 'GO_sites_genes_RAL_vs_FR_ZI_SAfr_Fst_autosomes.csv', 'GO_sites_genes_RAL_vs_FR_ZI_SAfr_Fst_ChrX.csv', 'GO_sites_genes_ZS_RAL_ZI_Fst_ChrX.csv']

#run in parallel
def run_in_parallel():
    pool = Pool(processes=18)
    pool.map(main_function, csv_list)


if __name__ == '__main__':
    run_in_parallel()




