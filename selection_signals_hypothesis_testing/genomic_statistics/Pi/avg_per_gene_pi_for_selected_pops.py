import pandas as pd


def find_gene_lengths(gene_map):

    #loading gene data
    genes = pd.read_csv(gene_map, keep_default_na=False)

    #finding lengths
    genes["length"] = genes["END"] - genes["START"] + 1

    #extracting lists and making dictionary
    fbgns = genes["FBgn"].to_list()
    lens = genes["length"].to_list()
    length_dict = dict(zip(fbgns, lens))

    return length_dict



def Pi_per_gene(csvfile, GeneMap):

    #loading data
    df = pd.read_csv(csvfile, keep_default_na=False)
    gene_lengths = find_gene_lengths(GeneMap)
    
    #comparison columns
    pop_cols = ["ZS_PI", "ZH_PI", "ZI_PI", "FR_PI", "RAL_PI"]

    #removing sites not mapped to a gene
    dfgenes = df[df['FBgn'] != "None"]

    #expanding sites with mutliple genes
    dfgenes.loc[:, 'FBgn'] = dfgenes['FBgn'].str.split('; ')
    dfgenes.loc[:, 'gene_symbol'] = dfgenes['gene_symbol'].str.split('; ')
    dfgenes_expand = dfgenes.explode(['FBgn', 'gene_symbol'])

    #averaging pi per gene
    dfgenes_expand = dfgenes_expand.drop(columns=['CHROM', 'POS'])
    pi_per_gene = dfgenes_expand.groupby(['FBgn', 'gene_symbol']).sum().reset_index()
    
    pi_per_gene["gene_length"] = pi_per_gene["FBgn"].map(gene_lengths)
    pi_per_gene = pi_per_gene.dropna(subset=['gene_length'])
    for col in pop_cols:
        pi_per_gene[col] = pi_per_gene[col] / pi_per_gene['gene_length']
    pi_per_gene.drop(columns=['gene_length'], inplace=True)


    return pi_per_gene







pipergene = Pi_per_gene("selected_pops_sites_pi_mapped2genes.csv", "Dmel_r6_protein_coding_genes.csv")

pipergene.to_csv("selected_pops_per_gene_pi.csv", index=False)


