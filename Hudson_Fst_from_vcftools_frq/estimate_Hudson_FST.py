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




def main_func(focal_pop):

    data = pd.read_csv("population_allele_frequencies.csv")

    pop_1 = focal_pop

    cos_pops = ["WCAfr", "EAfr", "Zam", "OOA-OW", "OOA-NW", "SAfr"]
    newcols = []


    for c in cos_pops:
        pop_2 = c

        fstcol, data = Hud_FST(data, pop_1, pop_2)
        newcols.append(fstcol)


    data = data[['CHROM', 'POS'] + newcols]

    data.to_csv(f"{pop_1}_fst.csv", index=False)




main_func("Zim")
main_func("ZS")
main_func("ZH")
main_func("ZC")
main_func("ZR")