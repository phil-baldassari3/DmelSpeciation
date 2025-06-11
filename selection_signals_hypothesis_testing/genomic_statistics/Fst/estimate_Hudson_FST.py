import pandas as pd

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

    data = pd.read_csv("selected_pops_allele_frequencies.csv")

    newcols = []

    for c in comparison_list:
        pop_1 = c[0]
        pop_2 = c[1]

        fstcol, data = Hud_FST(data, pop_1, pop_2)
        newcols.append(fstcol)


    data = data[['CHROM', 'POS'] + newcols]

    data.to_csv("selected_pops_fst.csv", index=False)





pop_comparisons = [("ZS", "RAL"), ("ZS", "FR"), ("ZS", "ZI"), ("ZH", "RAL"), ("ZH", "FR"), ("ZH", "ZI"), ("RAL", "FR"), ("FR", "ZI"), ("ZI", "RAL")]

main_func(pop_comparisons)