import locitableTK as lttk
from multiprocessing import Pool


def map_fst_csv2genes(csvfile):

    #load data
    fsttable = lttk.loci_table(csvfile, locitype="windowed", win_size=5000)

    #map genes
    fsttable.map2genes()

    #save new csv
    name = csvfile.replace(".csv", "_genes")
    fsttable.save(name)


#run in parallel
file_list = ["windowed_mean_5kbp_FST_averages_Autosomes.csv", "windowed_max_5kbp_FST_averages_ChromX.csv", "windowed_max_5kbp_FST_averages_Autosomes.csv", "windowed_mean_5kbp_FST_averages_ChromX.csv"]


def run_in_parallel():
    pool = Pool(processes=4)
    pool.map(map_fst_csv2genes, file_list)


if __name__ == '__main__':
    run_in_parallel()