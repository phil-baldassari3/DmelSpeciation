import os, sys
import pandas as pd
import numpy as np



####PARAMETERS TO SET##################################################
winsize = 10000
pct1 = 95
pct2 = 75
csv_list = []
#######################################################################

#functions
def find_idx(fst_ls, threshold):
    """
    Function takes in a list of Fst values and returns a list of the indeces of the values if they were greater than the given threshold.
    Returns: list of index values
    """

    idxls = [idx for idx, i in enumerate(fst_ls) if i > threshold]

    return idxls



def make_check_ls(idx_ls):
    """
    Function takes the idx list denoting positions from the Fst list that passed the threshold, and groups the positions into potential islands to be checked.
    Returns: list of 2-item lists
    """

    if len(idx_ls) == 1:
        out = [[idx_ls[0], idx_ls[0]]]

    else:
        out = [list(idx_ls[idx:idx+2]) for idx in range(len(idx_ls)-1)]

    return out



def convert2genomic(idx, win_size):
    """
    Takes an index integer and converts it to genomic position using the window size
    Returns: int
    """

    pos = (idx * win_size) + 1

    return pos



def island_idx2genomic(isidxls):
    """
    Function takes in the islands idx list and using the convert2genomic() function, converts the entire list of lists.
    It updates the list regardless of whether you assign output of the function to a variable or not.
    Returns: updated list of lists
    """

    for i in isidxls:
        for idx, j in enumerate(i):
            i[idx] = convert2genomic(j, winsize)

    return isidxls



def islands(checkls, fstls, threshold):
    """
    Function takes in the list of islands to check and the fst list which the check list came from.
    It then outputs a list of lists of islands that met the criteria given the threshold.
    Returns: List of lists
    """

    #first round island list
    islandls = []

    for i in checkls:
        if all(x > threshold for x in fstls[i[0]: i[1]]):
            islandls.append([i[0], i[1]+1])
        else:
            islandls.append([i[0], i[0]+1])

    #adding buffer at end
    islandls.append(-1)

    #combining adjacent islands (has to be repeated)
    for idx in range(len(islandls)):
        if islandls[idx] == -1 or islandls[idx+1] == -1:
            break
        elif islandls[idx] == []:
            continue
        else:
            next = 1
            while True:
                if islandls[idx+next] != -1 and islandls[idx][1] > islandls[idx+next][0]:
                    islandls[idx][1] = islandls[idx+next][1]
                    islandls[idx+next] = []
                    next += 1
                else:
                    break
            continue



    #getting rid of -1 and empty islands
    islandls.remove(-1)
    islandls = [i for i in islandls if i != []]

    return islandls


def df2megalist(df):
    """
    function takes a master windowed Fst dataframe and makes a list of every Fst value to be used for a distribution
    Returns: list
    """

    #dropping position columns
    df = df.drop(['Chrom', 'window_start', 'window_end'], axis=1)

    #megalist 
    megals = []

    #iterating over columns and outputting to mega list
    for i in df.columns:
       ls = df[i].to_list()
       megals += ls

    return megals



def df2dict(df):
    """
    Function takes the Fst dataframe and makes a dictionary of lists for each column. This dictionary is used to find the islands for each comparison
    Returns: disctionary of lists
    """

    #dropping position columns
    df = df.drop(['Chrom', 'window_start', 'window_end'], axis=1)

    #megalist 
    dictionary = {}

    #iterating over columns and outputting to mega list
    for i in df.columns:
       ls = df[i].to_list()
       dictionary.update({i: ls})

    return dictionary


    
def find_islands_from_df(df, thresh1, thresh2):
    """
    Function takes in a dataframe of pairwise Fst comparisons and finds islands for each pairwise comparison.
    thresh1 is the percentile threshold for selecting windows. thresh2 is the percentile threshold used to connect windowes.
    Return: dictionary of comparisons keys and list of island values
    """

    #finding thresholds
    megals = df2megalist(df)
    t1 = np.percentile(megals, thresh1)
    t2 = np.percentile(megals, thresh2)
    #print(t1,t2)

    #making comaprison dict
    comp_dict = df2dict(df)

    #finding islands for each comparison
    island_dict = {}

    for comparison in comp_dict.keys():
        idxes = find_idx(comp_dict[comparison], t1)
        check = make_check_ls(idxes)
        isl = islands(check, comp_dict[comparison], t2)
        isl = island_idx2genomic(isl)

        #update island dict
        island_dict.update({comparison: isl})

    return island_dict, t1, t2





def main(infile):
    """Main function. Takes in an input file and uses the rest of the functions to output csv files."""

    #opening windowed fst csv file (column for each comparison)
    df = pd.read_csv(infile)

    #finding islands
    islands_dict, threshold1, threshold2 = find_islands_from_df(df, pct1, pct2)

    #converting to df
    islands_df = pd.DataFrame.from_dict(islands_dict, orient='index')
    islands_df = islands_df.transpose()

    #saving file
    islands_df.to_csv(infile.replace('.csv', '_') + "islands.csv", index=False)

    #saving thresholds
    with open(infile.replace('.csv', '_') + "percentiles.txt", "w") as threshold_file:
        threshold_file.write("{pct}th Percentile: {fst}".format(pct=str(pct1), fst=str(threshold1)) + "\n{pct}th Percentile: {fst}".format(pct=str(pct2), fst=str(threshold2)))






main('win10kbp_ChrX_allcomparisons.csv')
main('win10kbp_Chr2L_allcomparisons.csv')
main('win10kbp_Chr2R_allcomparisons.csv')
main('win10kbp_Chr3L_allcomparisons.csv')
main('win10kbp_Chr3R_allcomparisons.csv')



