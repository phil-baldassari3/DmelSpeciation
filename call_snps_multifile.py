#3rd iteration of the python based SNP calling tool.  Goal is to make an alignment dataframe and filter out the monomorphic sites.

#importing modules
import pandas as pd
import numpy as np

#setting directory
directory = "/Users/philipbaldassari/Desktop/zim-cos_Ch3R"


########################################################################################################
print('reading zim fasta alignment file into python dictionary...')

#opening  zim fasta and converting to dictionary
with open("zim_3R.fa") as file_zim:
    fasta_zim = {line.strip(">\n"):next(file_zim).rstrip() for line in file_zim}
	
del file_zim

#converting values in zim dictionary to list
for k, v in fasta_zim.items():
	fasta_zim[k] = list(v)

print('making dataframe from zim dictionary...')

#making dataframe from dictionary
aln = pd.DataFrame.from_dict(fasta_zim)

del fasta_zim

########################################################################################################
print('reading ZI fasta alignment file into python dictionary...')

#opening ZI fasta and converting to dictionary
with open("ZI_3R.fa") as file_ZI:
    fasta_ZI = {line.strip(">\n"):next(file_ZI).rstrip() for line in file_ZI}
	
del file_ZI

#converting values in zim dictionary to list
for k, v in fasta_ZI.items():
	fasta_ZI[k] = list(v)
 
print('making dataframe from ZI dictionary...')

#adding ZI to the dataframe
ZI = pd.DataFrame.from_dict(fasta_ZI)

del fasta_ZI

#merging ZI to aln dataframe by index
aln = pd.merge(aln, ZI, left_index=True, right_index=True)

del ZI
########################################################################################################
print('reading SAfr fasta alignment file into python dictionary...')

#opening  zim fasta and converting to dictionary
with open("SAfr_3R.fa") as file_SAfr:
    fasta_SAfr = {line.strip(">\n"):next(file_SAfr).rstrip() for line in file_SAfr}
	
del file_SAfr

#converting values in zim dictionary to list
for k, v in fasta_SAfr.items():
	fasta_SAfr[k] = list(v)
 
print('making dataframe from SAfr dictionary...')

#making dataframe from dictionary
SAfr = pd.DataFrame.from_dict(fasta_SAfr)

del fasta_SAfr

#merging SAfr to aln dataframe by index
aln = pd.merge(aln, SAfr, left_index=True, right_index=True)

del SAfr
########################################################################################################
print('reading FR fasta alignment file into python dictionary...')

#opening  zim fasta and converting to dictionary
with open("FR_3R.fa") as file_FR:
    fasta_FR = {line.strip(">\n"):next(file_FR).rstrip() for line in file_FR}
	
del file_FR

#converting values in zim dictionary to list
for k, v in fasta_FR.items():
	fasta_FR[k] = list(v)
 
print('making dataframe from FR dictionary...')

#making dataframe from dictionary
FR = pd.DataFrame.from_dict(fasta_FR)

del fasta_FR

#merging SAfr to aln dataframe by index
aln = pd.merge(aln, FR, left_index=True, right_index=True)

del FR
########################################################################################################
print('reading RAL fasta alignment file into python dictionary...')

#opening  zim fasta and converting to dictionary
with open("RAL_3R.fa") as file_RAL:
    fasta_RAL = {line.strip(">\n"):next(file_RAL).rstrip() for line in file_RAL}
	
del file_RAL

#converting values in zim dictionary to list
for k, v in fasta_RAL.items():
	fasta_RAL[k] = list(v)
 
print('making dataframe from RAL dictionary...')

#making dataframe from dictionary
RAL = pd.DataFrame.from_dict(fasta_RAL)

del fasta_RAL

#merging RAL to aln dataframe by index
aln = pd.merge(aln, RAL, left_index=True, right_index=True)

del RAL

print('Done making preliminatry alignment dataframe...')
########################################################################################################
print('replacing missing nucleotides with NaN...')
#Replace N with NaN
aln = aln.replace('N', np.nan)

print('labeling loci...')
#Add Locus column
aln.insert(0, "Locus", aln.index + 1)

print('done...continuing')
print('counting alleles...')

#new dataframe for bp counts
counts = pd.DataFrame(columns = ["Acount", "Tcount", "Ccount", "Gcount"])

#iterating through aln dataframe and counting bps for counts dataframe by tuples
for row in aln.itertuples():
                
        Acount = row.count("A")
        Tcount = row.count("T")
        Ccount = row.count("C")
        Gcount = row.count("G")
        
        counts = counts.append({'Acount' : Acount, 'Tcount' : Tcount, 'Ccount' : Ccount, 'Gcount' : Gcount}, 
                ignore_index = True)

#merging counts to aln dataframe by index
aln = pd.merge(aln, counts, left_index=True, right_index=True)

del counts

print('done...continuing')
print('estimating allele frequecies at each site...')

#summing for total bp count
aln["bpcount"] = aln['Acount'] + aln['Tcount'] + aln['Ccount'] + aln['Gcount']

#getting rid of rows with zero bps
aln = aln[aln.bpcount !=0]

#allele frqeuncies per nucleotide
aln['Aprop'] = aln['Acount']/aln['bpcount']
aln['Tprop'] = aln['Tcount']/aln['bpcount']
aln['Cprop'] = aln['Ccount']/aln['bpcount']
aln['Gprop'] = aln['Gcount']/aln['bpcount']

#Major allele Frequency
aln['MajorAF'] = aln[["Aprop","Tprop","Cprop","Gprop"]].max(axis=1)

print('done...continuing')
print('getting rid of monomorphic sites...')

#getting rid of monomorphic sites
aln = aln[aln.MajorAF !=1]

print('done...finalizing')

#saving as csv
aln.to_csv('zim-cos_Chr3R_SNPs.csv', index=False)

print('Process complete! A file has been saved to this directory named "zim-cos_Chr3R_SNPs.csv"')
  
