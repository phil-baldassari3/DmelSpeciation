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
file_list = ["Zim_fst.csv", "ZS_fst.csv", "ZH_fst.csv", "ZC_fst.csv", "ZR_fst.csv"]

def run_in_parallel():
    pool = Pool(processes=5)
    pool.map(map_merged_fst_csv2genes, file_list)


if __name__ == '__main__':
    run_in_parallel()