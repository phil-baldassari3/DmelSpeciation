import locitableTK as lttk
from multiprocessing import Pool


def map_merged_fst_csv2genes(csvfile):

    #load data
    fsttable = lttk.loci_table(csvfile)

    #map genes
    fsttable.map2genes()

    #save new csv
    name = csvfile.replace(".csv", "_genes")
    fsttable.save(name)


#run in parallel
file_list = ["ZR_diff_SNPs_Autosomes.csv", "ZC_diff_SNPs_Autosomes.csv", 
             "ZH_diff_SNPs_Autosomes.csv", "ZS_diff_SNPs_Autosomes.csv", 
             "Zim_diff_SNPs_Autosomes.csv", "ZR_diff_SNPs_ChromX.csv", 
             "ZC_diff_SNPs_ChromX.csv", "ZH_diff_SNPs_ChromX.csv", 
             "ZS_diff_SNPs_ChromX.csv", "Zim_diff_SNPs_ChromX.csv"]


def run_in_parallel():
    pool = Pool(processes=10)
    pool.map(map_merged_fst_csv2genes, file_list)


if __name__ == '__main__':
    run_in_parallel()