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



def FST_per_gene(csvfile, GeneMap):

    #loading data
    df = pd.read_csv(csvfile, keep_default_na=False)
    gene_lengths = find_gene_lengths(GeneMap)
    
    #comparison columns
    comparison_cols = ["ZS.vs.RAL_FST", "ZS.vs.FR_FST", "ZS.vs.ZI_FST", "ZH.vs.RAL_FST", "ZH.vs.FR_FST", "ZH.vs.ZI_FST", "RAL.vs.FR_FST", "FR.vs.ZI_FST", "ZI.vs.RAL_FST"]

    #removing sites not mapped to a gene
    dfgenes = df[df['FBgn'] != "None"]

    #clipping Fst values to zero
    dfgenes.loc[:, comparison_cols] = dfgenes.loc[:, comparison_cols].clip(lower=0)

    #expanding sites with mutliple genes
    dfgenes.loc[:, 'FBgn'] = dfgenes['FBgn'].str.split('; ')
    dfgenes.loc[:, 'gene_symbol'] = dfgenes['gene_symbol'].str.split('; ')
    dfgenes_expand = dfgenes.explode(['FBgn', 'gene_symbol'])

    #averaging pi per gene
    dfgenes_expand = dfgenes_expand.drop(columns=['CHROM', 'POS'])
    fst_per_gene = dfgenes_expand.groupby(['FBgn', 'gene_symbol']).sum().reset_index()
    
    fst_per_gene["gene_length"] = fst_per_gene["FBgn"].map(gene_lengths)
    fst_per_gene = fst_per_gene.dropna(subset=['gene_length'])
    for col in comparison_cols:
        fst_per_gene[col] = fst_per_gene[col] / fst_per_gene['gene_length']
    fst_per_gene.drop(columns=['gene_length'], inplace=True)


    return fst_per_gene







fstpergene = FST_per_gene("selected_pops_persite_fst_mapped2genes.csv", "Dmel_r6_protein_coding_genes.csv")

fstpergene.to_csv("selected_pops_per_gene_fst_NormByGeneLength.csv", index=False)


