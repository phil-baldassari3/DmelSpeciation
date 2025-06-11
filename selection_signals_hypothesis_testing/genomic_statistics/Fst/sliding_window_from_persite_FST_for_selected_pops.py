import pandas as pd


def chromosome_splitter(df):
    """
    splits an autosome df into searate chromosome dataframes
    """

    df_2L = df[df["CHROM"] == "2L"]
    df_2R = df[df["CHROM"] == "2R"]
    df_3L = df[df["CHROM"] == "3L"]
    df_3R = df[df["CHROM"] == "3R"]
    df_X = df[df["CHROM"] == "X"]

    return df_2L, df_2R, df_3L, df_3R, df_X



def sliding_window(df, step, sumstat):
    """
    Takes df of per site Fst and returns a sliding window Fst. If this is an autosome dataframe you need to run it through autosome_splitter() first
    Can compute average per site FST per window or max per site FST per window based on the sumstat argument which can be either "mean" or "max"
    """

    #setting chrom
    chrom_ls = df['CHROM'].tolist()

    if 'X' in chrom_ls:
        chrom = 'X'
        chrom_len = 23542271
    elif '2L' in chrom_ls:
        chrom = '2L'
        chrom_len = 23513712
    elif '2R' in chrom_ls:
        chrom = '2R'
        chrom_len = 25286936
    elif '3L' in chrom_ls:
        chrom = '3L'
        chrom_len = 28110227
    elif '3R' in chrom_ls:
        chrom = '3R'
        chrom_len = 32079331
    
    
    #setting avg fst lists
    win_fst_ZSvsRAL = []
    win_fst_ZSvsFR = []
    win_fst_ZSvsZI = []
    win_fst_ZHvsRAL = []
    win_fst_ZHvsFR = []
    win_fst_ZHvsZI = []
    win_fst_RALvsFR = []
    win_fst_FRvsZI = []
    win_fst_ZIvsRAL = []

    #window
    win_start = []
    win_end = []


    for i in range(0, chrom_len, step):

        #window lists
        win_start.append(i+1)
        win_end.append(i+step)

        #chunking df
        chunk_df = df[(df["POS"] >= i+1) & (df["POS"] < i+step)]


        if len(chunk_df.index) == 0:
            win_fst_ZSvsRAL.append(0)
            win_fst_ZSvsFR.append(0)
            win_fst_ZSvsZI.append(0)
            win_fst_ZHvsRAL.append(0)
            win_fst_ZHvsFR.append(0)
            win_fst_ZHvsZI.append(0)
            win_fst_RALvsFR.append(0)
            win_fst_FRvsZI.append(0)
            win_fst_ZIvsRAL.append(0)

        else:

            if sumstat == "mean":
                winZSvsRAL = chunk_df.iloc[:, 2].mean()
                winZSvsFR = chunk_df.iloc[:, 3].mean()
                winZSvsZI = chunk_df.iloc[:, 4].mean()
                winZHvsRAL = chunk_df.iloc[:, 5].mean()
                winZHvsFR = chunk_df.iloc[:, 6].mean()
                winZHvsZI = chunk_df.iloc[:, 7].mean()
                winRALvsFR = chunk_df.iloc[:, 8].mean()
                winFRvsZI = chunk_df.iloc[:, 9].mean()
                winZIvsRAL = chunk_df.iloc[:, 10].mean()
            elif sumstat == "max":
                winZSvsRAL = chunk_df.iloc[:, 2].max()
                winZSvsFR = chunk_df.iloc[:, 3].max()
                winZSvsZI = chunk_df.iloc[:, 4].max()
                winZHvsRAL = chunk_df.iloc[:, 5].max()
                winZHvsFR = chunk_df.iloc[:, 6].max()
                winZHvsZI = chunk_df.iloc[:, 7].max()
                winRALvsFR = chunk_df.iloc[:, 8].max()
                winFRvsZI = chunk_df.iloc[:, 9].max()
                winZIvsRAL = chunk_df.iloc[:, 10].max()
            else:
                print("wrong sumstat argument")
                break


            win_fst_ZSvsRAL.append(winZSvsRAL)
            win_fst_ZSvsFR.append(winZSvsFR)
            win_fst_ZSvsZI.append(winZSvsZI)
            win_fst_ZHvsRAL.append(winZHvsRAL)
            win_fst_ZHvsFR.append(winZHvsFR)
            win_fst_ZHvsZI.append(winZHvsZI)
            win_fst_RALvsFR.append(winRALvsFR)
            win_fst_FRvsZI.append(winFRvsZI)
            win_fst_ZIvsRAL.append(winZIvsRAL)



        dictionary = {
            "START": win_start,
            "END": win_end,
            "Win_FST_ZS.vs.RAL" : win_fst_ZSvsRAL,
            "Win_FST_ZS.vs.FR" : win_fst_ZSvsFR,
            "Win_FST_ZS.vs.ZI" : win_fst_ZSvsZI,
            "Win_FST_ZH.vs.RAL" : win_fst_ZHvsRAL,
            "Win_FST_ZH.vs.FR" : win_fst_ZHvsFR,
            "Win_FST_ZH.vs.ZI" : win_fst_ZHvsZI,
            "Win_FST_RAL.vs.FR" : win_fst_RALvsFR,
            "Win_FST_FR.vs.ZI" : win_fst_FRvsZI,
            "Win_FST_ZI.vs.RAL" : win_fst_ZIvsRAL
        }

        windowed_df = pd.DataFrame.from_dict(dictionary)


    #adding chrom column
    windowed_df.insert(0, 'CHROM', chrom)

    return windowed_df




def main_function(infile, window, summary_statistic):
    """
    main function of the program. uses sliding_window() and autosome_splitter() to output annotated sliding window Fst files for each infile. THis function in run in parallel for each infile
    """

    persite_df = pd.read_csv(infile)
    comparison_cols = ["ZS.vs.RAL_FST", "ZS.vs.FR_FST", "ZS.vs.ZI_FST", "ZH.vs.RAL_FST", "ZH.vs.FR_FST", "ZH.vs.ZI_FST", "RAL.vs.FR_FST", "FR.vs.ZI_FST", "ZI.vs.RAL_FST"]
    persite_df.loc[:, comparison_cols] = persite_df.loc[:, comparison_cols].clip(lower=0)


    persite_df_2L, persite_df_2R, persite_df_3L, persite_df_3R, persite_df_X = chromosome_splitter(persite_df)

    win_df_2L = sliding_window(persite_df_2L, window, summary_statistic)
    win_df_2R = sliding_window(persite_df_2R, window, summary_statistic)
    win_df_3L = sliding_window(persite_df_3L, window, summary_statistic)
    win_df_3R = sliding_window(persite_df_3R, window, summary_statistic)
    win_df_X = sliding_window(persite_df_X, window, summary_statistic)

    win_df = pd.concat([win_df_2L, win_df_2R, win_df_3L, win_df_3R, win_df_X])
    win_df.to_csv('windowed_{st}_'.format(st=summary_statistic) + str(int(window/1000)) + 'kbp_{file}'.format(file=infile), index=False)






main_function("selected_pops_persite_fst.csv", 10000, "mean")


