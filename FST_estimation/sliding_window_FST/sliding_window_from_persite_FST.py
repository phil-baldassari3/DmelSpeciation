import pandas as pd


def autosome_splitter(df):
    """
    splits an autosome df into searate chromosome dataframes
    """

    df_2L = df[df["CHROM"] == "2L"]
    df_2R = df[df["CHROM"] == "2R"]
    df_3L = df[df["CHROM"] == "3L"]
    df_3R = df[df["CHROM"] == "3R"]

    return df_2L, df_2R, df_3L, df_3R



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
    win_fst_Zim = []
    win_fst_ZS = []
    win_fst_ZH = []
    win_fst_ZC = []
    win_fst_ZR = []
    win_fst_Cos = []

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
            win_fst_Zim.append(0)
            win_fst_ZS.append(0)
            win_fst_ZH.append(0)
            win_fst_ZC.append(0)
            win_fst_ZR.append(0)
            win_fst_Cos.append(0)

        else:

            if sumstat == "mean":
                winZim = chunk_df.iloc[:, 2].mean()
                winZS = chunk_df.iloc[:, 3].mean()
                winZH = chunk_df.iloc[:, 4].mean()
                winZC = chunk_df.iloc[:, 5].mean()
                winZR = chunk_df.iloc[:, 6].mean()
                winCos = chunk_df.iloc[:, 7].mean()
            elif sumstat == "max":
                winZim = chunk_df.iloc[:, 2].max()
                winZS = chunk_df.iloc[:, 3].max()
                winZH = chunk_df.iloc[:, 4].max()
                winZC = chunk_df.iloc[:, 5].max()
                winZR = chunk_df.iloc[:, 6].max()
                winCos = chunk_df.iloc[:, 7].max()
            else:
                print("wrong sumstat argument")
                break


            win_fst_Zim.append(winZim)
            win_fst_ZS.append(winZS)
            win_fst_ZH.append(winZH)
            win_fst_ZC.append(winZC)
            win_fst_ZR.append(winZR)
            win_fst_Cos.append(winCos)



        dictionary = {
            "window_start": win_start,
            "window_end": win_end,
            "Win_FST_Zim": win_fst_Zim,
            "Win_FST_ZS": win_fst_ZS,
            "Win_FST_ZH": win_fst_ZH,
            "Win_FST_ZC": win_fst_ZC,
            "Win_FST_ZR": win_fst_ZR,
            "Win_FST_Cos": win_fst_Cos
        }

        windowed_df = pd.DataFrame.from_dict(dictionary)


    #adding chrom column
    windowed_df.insert(0, 'Chrom', chrom)

    return windowed_df




def main_function(infile, window, summary_statistic):
    """
    main function of the program. uses sliding_window() and autosome_splitter() to output annotated sliding window Fst files for each infile. THis function in run in parallel for each infile
    """

    persite_df = pd.read_csv(infile)

    if "ChromX" in infile:

        win_df = sliding_window(persite_df, window, summary_statistic)
        win_df.to_csv('windowed_{st}_'.format(st=summary_statistic) + str(int(window/1000)) + 'kbp_{file}'.format(file=infile), index=False)


    else:

        persite_df_2L, persite_df_2R, persite_df_3L, persite_df_3R = autosome_splitter(persite_df)

        win_df_2L = sliding_window(persite_df_2L, window, summary_statistic)
        win_df_2R = sliding_window(persite_df_2R, window, summary_statistic)
        win_df_3L = sliding_window(persite_df_3L, window, summary_statistic)
        win_df_3R = sliding_window(persite_df_3R, window, summary_statistic)

        win_df = pd.concat([win_df_2L, win_df_2R, win_df_3L, win_df_3R])
        win_df.to_csv('windowed_{st}_'.format(st=summary_statistic) + str(int(window/1000)) + 'kbp_{file}'.format(file=infile), index=False)






main_function("FST_averages_Autosomes.csv", 5000, "mean")
main_function("FST_averages_ChromX.csv", 5000, "mean")
main_function("FST_averages_Autosomes.csv", 5000, "max")
main_function("FST_averages_ChromX.csv", 5000, "max")