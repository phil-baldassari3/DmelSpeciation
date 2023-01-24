#importing modules
import os,sys
import pandas as pd
import numpy as np
from statistics import mean
from statistics import stdev
import scipy.stats
import matplotlib.pyplot as plt

from multiprocessing import Pool

from annotate_ztest_toolkit import loci_table




#setting directory
directory = "/Users/philipbaldassari/Desktop/zim-cos_gene_GO_ann/10kbp_windowed_GO_ann"


#get list of txt list files for parallel pool mapping
file_list = []

for file in os.listdir(directory):
    if file.endswith('.csv'):
        file_list.append(file)
    else:
        continue

print(file_list)
print(len(file_list))
print('\n\n\n')




#enrichment funciton using toolkit
def win_enrich_with_toolkit(infile):

    print('starting on' + infile)

    #loci table instance (no need to map and annotate with these files)
    fstwins1 = loci_table("10kbp_windowed_GO_ann/" + infile, locitype="windowed", win_size=10000)
    fstwins5 = loci_table("10kbp_windowed_GO_ann/" + infile, locitype="windowed", win_size=10000)
    #fstwins10 = loci_table("10kbp_windowed_GO_ann/" + infile, locitype="windowed", win_size=10000)

    #performing the enrichment test
    print("starting neurogenesis enrichment, threshold 1%" + str(infile))
    fstwins1.enrichment_test(10000, 'Avg_Fst', 1, 'Neurogenesis', 'GO_gene_lists/neurogenesis_genes.txt')
    print("starting Mating_behavior enrichment, threshold 1%" + str(infile))
    fstwins1.enrichment_test(10000, 'Avg_Fst', 1, 'Mating_behavior', 'GO_gene_lists/mating_behavior_genes.txt')
    print("starting Male_mating_behavior enrichment, threshold 1%" + str(infile))
    fstwins1.enrichment_test(10000, 'Avg_Fst', 1, 'Male_mating_behavior', 'GO_gene_lists/male_mating_genes.txt')
    print("starting Female_mating_behavior enrichment, threshold 1%" + str(infile))
    fstwins1.enrichment_test(10000, 'Avg_Fst', 1, 'Female_mating_behavior', 'GO_gene_lists/female_mating_genes.txt')


    print("starting neurogenesis enrichment, threshold 5%" + str(infile))
    fstwins5.enrichment_test(10000, 'Avg_Fst', 5, 'Neurogenesis', 'GO_gene_lists/neurogenesis_genes.txt')
    print("starting Mating_behavior enrichment, threshold 5%" + str(infile))
    fstwins5.enrichment_test(10000, 'Avg_Fst', 5, 'Mating_behavior', 'GO_gene_lists/mating_behavior_genes.txt')
    print("starting Male_mating_behavior enrichment, threshold 5%" + str(infile))
    fstwins5.enrichment_test(10000, 'Avg_Fst', 5, 'Male_mating_behavior', 'GO_gene_lists/male_mating_genes.txt')
    print("starting Female_mating_behavior enrichment, threshold 5%" + str(infile))
    fstwins5.enrichment_test(10000, 'Avg_Fst', 5, 'Female_mating_behavior', 'GO_gene_lists/female_mating_genes.txt')

    """ 
    print("starting neurogenesis enrichment, threshold 10%" + str(infile))
    fstwins10.enrichment_test(10000, 'Avg_Fst', 10, 'Neurogenesis', 'GO_gene_lists/neurogenesis_genes.txt')
    print("starting Mating_behavior enrichment, threshold 10%" + str(infile))
    fstwins10.enrichment_test(10000, 'Avg_Fst', 10, 'Mating_behavior', 'GO_gene_lists/mating_behavior_genes.txt')
    print("starting Male_mating_behavior enrichment, threshold 10%" + str(infile))
    fstwins10.enrichment_test(10000, 'Avg_Fst', 10, 'Male_mating_behavior', 'GO_gene_lists/male_mating_genes.txt')
    print("starting Female_mating_behavior enrichment, threshold 10%" + str(infile))
    fstwins10.enrichment_test(10000, 'Avg_Fst', 10, 'Female_mating_behavior', 'GO_gene_lists/female_mating_genes.txt')
    """

    #printing Z test resuts to file
    original_stdout = sys.stdout


    with open("windowed_enrichment_analysis/" + "results_" + str(infile) + "fstwin_1percent.txt", 'w') as file:
        sys.stdout = file 
        print(fstwins1)
        sys.stdout = original_stdout

    with open("windowed_enrichment_analysis/" + "results_" + str(infile) + "fstwin_5percent.txt", 'w') as file:
        sys.stdout = file 
        print(fstwins5)
        sys.stdout = original_stdout

    """ with open("windowed_enrichment_analysis/" + "results_" + str(infile) + "fstwin_10percent.txt", 'w') as file:
        sys.stdout = file 
        print(fstwins10)
        sys.stdout = original_stdout """

    

    #distribution plots
    namedist = "windowed_enrichment_analysis/dist_Fst" + str(infile)
    fstwins1.plot_dist_param("Avg_Fst", namedist)


    #z test plots
    fstwins1.plot_enrichment_test("windowed_enrichment_analysis/" + str(infile) + "neurogenesis_win_top1", "Neurogenesis")
    fstwins1.plot_enrichment_test("windowed_enrichment_analysis/" + str(infile) + "Mating_behavior_win_top1", "Mating_behavior")
    fstwins1.plot_enrichment_test("windowed_enrichment_analysis/" + str(infile) + "Male_mating_behavior_win_top1", "Male_mating_behavior")
    fstwins1.plot_enrichment_test("windowed_enrichment_analysis/" + str(infile) + "Female_mating_behavior_win_top1", "Female_mating_behavior")


    fstwins5.plot_enrichment_test("windowed_enrichment_analysis/" + str(infile) + "neurogenesis_win_top5", "Neurogenesis")
    fstwins5.plot_enrichment_test("windowed_enrichment_analysis/" + str(infile) + "Mating_behavior_win_top5", "Mating_behavior")
    fstwins5.plot_enrichment_test("windowed_enrichment_analysis/" + str(infile) + "Male_mating_behavior_win_top5", "Male_mating_behavior")
    fstwins5.plot_enrichment_test("windowed_enrichment_analysis/" + str(infile) + "Female_mating_behavior_win_top5", "Female_mating_behavior")

    """ 
    fstwins10.plot_enrichment_test("windowed_enrichment_analysis/" + str(infile) + "neurogenesis_win_top10", "Neurogenesis")
    fstwins10.plot_enrichment_test("windowed_enrichment_analysis/" + str(infile) + "Mating_behavior_win_top10", "Mating_behavior")
    fstwins10.plot_enrichment_test("windowed_enrichment_analysis/" + str(infile) + "Male_mating_behavior_win_top10", "Male_mating_behavior")
    fstwins10.plot_enrichment_test("windowed_enrichment_analysis/" + str(infile) + "Female_mating_behavior_win_top10", "Female_mating_behavior")
    """




#for parallel mapping

#run in parallel
def run_in_parallel():
    pool = Pool(processes=18)
    pool.map(win_enrich_with_toolkit, file_list)


if __name__ == '__main__':
    run_in_parallel()

'''





test = loci_table('GO_sites_genes_ZS_RAL_ZI_Fst_autosomes.csv', locitype="windowed", win_size=10000)

test.save("test")

'''