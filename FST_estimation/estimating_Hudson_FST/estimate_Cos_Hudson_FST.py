import pandas as pd
from itertools import combinations

#functions
def est_FST(row, n1, n2, p1, p2):
    
    try:
        fst = (((row[p1] - row[p2])**2) - ((row[p1]*(1-row[p1]))/(row[n1]-1)) - ((row[p2]*(1-row[p2]))/(row[n2]-1))) / ((row[p1]*(1-row[p2])) + (row[p2]*(1-row[p1])))
    except ZeroDivisionError:
        fst = 0

    return fst


def Hud_FST(df, pop1, pop2):

    newcol = f"{pop1}.vs.{pop2}_FST"

    df[newcol] = df.apply(lambda row: est_FST(row, f'{pop1}_N', f'{pop2}_N', f'{pop1}_maf', f'{pop2}_maf'), axis=1)

    return newcol, df




def main_func(comparison_list):

    data = pd.read_csv("population_allele_frequencies.csv")

    newcols = []

    for c in comparison_list:
        pop_1 = c[0]
        pop_2 = c[1]

        fstcol, data = Hud_FST(data, pop_1, pop_2)
        newcols.append(fstcol)


    data = data[['CHROM', 'POS'] + newcols]

    data.to_csv("Cos_fst.csv", index=False)





cos_pops = ["WCAfr", "EAfr", "Zam", "OOA-OW", "OOA-NW", "SAfr"]
cos_combinations = list(combinations(cos_pops, 2))


main_func(cos_combinations)