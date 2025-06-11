import pandas as pd
import numpy as np

#set window size
window_size = 10000




#gene locations from ensembl biomart
genes = pd.read_csv("Dmel_r6_protein_coding_genes.csv", keep_default_na=False)

#windowed FST data
winFST = pd.read_csv("windowed_mean_10kbp_selected_pops_fst.csv")



##### preparing gene df for merge #####

#rounding START down to nearest window_size + 1 and rounding END up to nearest window_size
genes['START'] = (genes['START'] // window_size) * window_size + 1
genes['END'] = ((genes['END'] + (window_size - 1)) // window_size) * window_size

#split intervals into windows
genes['winstarts'] = genes.apply(lambda row: list(range(row['START'], row['END'], window_size)), axis=1)
genes = genes.explode('winstarts', ignore_index=True)
genes['START'] = genes['winstarts']
genes['END'] = genes['START'] + window_size - 1
genes = genes.drop(columns=['winstarts'])
genes = genes.groupby(["CHROM", "START", "END"], as_index=False).agg({"FBgn": lambda x: "; ".join(x), "gene_symbol": lambda x: "; ".join(x)})


##### merging with windowed FST data #####

#merge on CHROM START END
merged = winFST.merge(genes, on=["CHROM", "START", "END"], how="left")
merged["FBgn"] = merged["FBgn"].fillna("None")
merged["gene_symbol"] = merged["gene_symbol"].fillna("None")


#saving csv
merged.to_csv("windowed_mean_10kbp_selected_pops_fst_mapped2genes.csv", index=False)