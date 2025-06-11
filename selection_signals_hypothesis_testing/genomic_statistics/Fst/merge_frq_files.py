import pandas as pd
import os


def load_allele_freq_data(frqfile):


    #Function to extract the minimum allele frequency
    def extract_min_maf(allele_str):
        freq = float(allele_str.split(':')[1])
        return freq
    


    popname = frqfile.replace(".frq", "")

    #load data
    df = pd.read_csv(f"{popname}.frq", sep='\t')

    #Drop the 'N_ALLELES' column
    df = df.drop(columns=['N_ALLELES', f'MAF_{popname}'])

    #Apply the function to maf column
    df[f'{popname}_maf'] = df[f'{popname}_maf'].apply(extract_min_maf)

    return df




#df list
dfs = []

for file in os.listdir():
    if file.endswith(".frq"):
        
        df = load_allele_freq_data(file)
        dfs.append(df)


#starting df
merged_df = dfs[0]

#looping and merging remaining dfs
for df in dfs[1:]:
    merged_df = pd.merge(merged_df, df, on=["CHROM", "POS"], how="inner")


merged_df.to_csv("selected_pops_allele_frequencies.csv", index=False)