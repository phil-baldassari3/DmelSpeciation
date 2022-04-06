#takes vcf_Chr_samples-sites.csv files and makes .txt filtration files to be used with vcftools --positions

#importing modules
import numpy as np
import pandas as pd
from multiprocessing import Process

def filtrationX():
    print("finding quality sites and filtering maf 5% for ChrX")

    #open dataframe
    df = pd.read_csv("vcf_ChrX_samples-sites.csv")

    #estimating maf and sample coverage
    df["Total Allele"] = df["Major Allele Count"] + df["Minor Allele Count"]

    df["maf"] = df["Minor Allele Count"]/df["Total Allele"]

    df["Sample_Coverage"] = df["Total Allele"]/(df["Total Allele"] + df["Missing Sample Count"])

    #remvoing sites with >45% coverage
    df = df[df.Sample_Coverage >=0.45]

    #maf 5%
    df = df[df.maf >=0.05]

    #making filtration file
    df = df.drop(['Major Allele Count', 'Minor Allele Count', 'Missing Sample Count', 'Total Allele', 'maf', 'Sample_Coverage'], axis=1)

    df.to_csv('vcf_ChrX_sites-filtration.txt', sep="\t", index=False)

    print("saved ChrX vcf filtration file")


def filtration2L():
    print("finding quality sites and filtering maf 5% for Chr2L")

    #open dataframe
    df = pd.read_csv("vcf_Chr2L_samples-sites.csv")

    #estimating maf and sample coverage
    df["Total Allele"] = df["Major Allele Count"] + df["Minor Allele Count"]

    df["maf"] = df["Minor Allele Count"]/df["Total Allele"]

    df["Sample_Coverage"] = df["Total Allele"]/(df["Total Allele"] + df["Missing Sample Count"])

    #remvoing sites with >45% coverage
    df = df[df.Sample_Coverage >=0.45]

    #maf 5%
    df = df[df.maf >=0.05]

    #making filtration file
    df = df.drop(['Major Allele Count', 'Minor Allele Count', 'Missing Sample Count', 'Total Allele', 'maf', 'Sample_Coverage'], axis=1)

    df.to_csv('vcf_Chr2L_sites-filtration.txt', sep="\t", index=False)

    print("saved Chr2L vcf filtration file")


def filtration2R():
    print("finding quality sites and filtering maf 5% for Chr2R")

    #open dataframe
    df = pd.read_csv("vcf_Chr2R_samples-sites.csv")

    #estimating maf and sample coverage
    df["Total Allele"] = df["Major Allele Count"] + df["Minor Allele Count"]

    df["maf"] = df["Minor Allele Count"]/df["Total Allele"]

    df["Sample_Coverage"] = df["Total Allele"]/(df["Total Allele"] + df["Missing Sample Count"])

    #remvoing sites with >45% coverage
    df = df[df.Sample_Coverage >=0.45]

    #maf 5%
    df = df[df.maf >=0.05]

    #making filtration file
    df = df.drop(['Major Allele Count', 'Minor Allele Count', 'Missing Sample Count', 'Total Allele', 'maf', 'Sample_Coverage'], axis=1)

    df.to_csv('vcf_Chr2R_sites-filtration.txt', sep="\t", index=False)

    print("saved Chr2R vcf filtration file")


def filtration3L():
    print("finding quality sites and filtering maf 5% for Chr3L")

    #open dataframe
    df = pd.read_csv("vcf_Chr3L_samples-sites.csv")

    #estimating maf and sample coverage
    df["Total Allele"] = df["Major Allele Count"] + df["Minor Allele Count"]

    df["maf"] = df["Minor Allele Count"]/df["Total Allele"]

    df["Sample_Coverage"] = df["Total Allele"]/(df["Total Allele"] + df["Missing Sample Count"])

    #remvoing sites with >45% coverage
    df = df[df.Sample_Coverage >=0.45]

    #maf 5%
    df = df[df.maf >=0.05]

    #making filtration file
    df = df.drop(['Major Allele Count', 'Minor Allele Count', 'Missing Sample Count', 'Total Allele', 'maf', 'Sample_Coverage'], axis=1)

    df.to_csv('vcf_Chr3L_sites-filtration.txt', sep="\t", index=False)

    print("saved Chr3L vcf filtration file")


def filtration3R():
    print("finding quality sites and filtering maf 5% for Chr3R")

    #open dataframe
    df = pd.read_csv("vcf_Chr3R_samples-sites.csv")

    #estimating maf and sample coverage
    df["Total Allele"] = df["Major Allele Count"] + df["Minor Allele Count"]

    df["maf"] = df["Minor Allele Count"]/df["Total Allele"]

    df["Sample_Coverage"] = df["Total Allele"]/(df["Total Allele"] + df["Missing Sample Count"])

    #remvoing sites with >45% coverage
    df = df[df.Sample_Coverage >=0.45]

    #maf 5%
    df = df[df.maf >=0.05]

    #making filtration file
    df = df.drop(['Major Allele Count', 'Minor Allele Count', 'Missing Sample Count', 'Total Allele', 'maf', 'Sample_Coverage'], axis=1)

    df.to_csv('vcf_Chr3R_sites-filtration.txt', sep="\t", index=False)

    print("saved Chr3R vcf filtration file")


#Running concurrently

if __name__ == '__main__':
	pX = Process(target=filtrationX)
	pX.start()
	
	p2L = Process(target=filtration2L)
	p2L.start()
	
	p2R = Process(target=filtration2R)
	p2R.start()
	
	p3L = Process(target=filtration3L)
	p3L.start()
	
	p3R = Process(target=filtration3R)
	p3R.start()
	
	pX.join()
	p2L.join()
	p2R.join()
	p3L.join()
	p3R.join()




