#importing modules
import pandas as pd
import numpy as np
from Bio import SeqIO


#setting directory
directory = "/Users/philipbaldassari/desktop/fst_test/"

#getting length of sequence
seq = open('dummy.seq')
seq_data = seq.read()
seqlen = len(seq_data)
seq.close()

#set window size
win = 7

#counter for loop, used for loci data
counter = 0

#setting window with range function
for i in range(0, seqlen, win):
    j = seqlen if i+win>seqlen else i+win

    #dataframes from fasta using Seq.IO
    records = SeqIO.parse("aln_test.fa", 'fasta')
    seqList = []
    for record in records:
        desp = record.description
        seq = list(record.seq._data.upper()[i:j]) #indeces from range window
        seqList.append([desp] + seq)
        aln = pd.DataFrame(seqList)
        aln = aln.T
        aln.columns = aln.iloc[0]
        aln = aln.drop(aln.index[0])
        aln = aln.reset_index(drop=True)
    #print(aln)


    #Replace N with NaN
    aln = aln.replace('N', np.nan)

    #Add Locus column
    aln.insert(0, "Locus", aln.index + (win*counter) + 1)

    #changing the counter
    counter += 1

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


    #getting rid of monomorphic sites
    aln = aln[aln.MajorAF !=1]


    #saving as csv
    aln.to_csv('SNPs_called_{num}.csv'.format(num=counter), index=False)


print('Done!')  

