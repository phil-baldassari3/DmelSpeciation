import pandas as pd



def FST_per_gene_across_comparisons(csvfile, focal_pop, sumstat):
    """
    sumstat can be set to either "mean" or "max"
    """

    #loading data
    df = pd.read_csv(csvfile)


    
    #optional filtering
    print(df)
    df = df[df['CHROM'] == "X"]
    print(df)
    
    



    #comparison columns
    comparison_cols = [f'{focal_pop}.vs.WCAfr_FST',
                    f'{focal_pop}.vs.SAfr_FST',
                    f'{focal_pop}.vs.OOA-NW_FST',
                    f'{focal_pop}.vs.EAfr_FST',
                    f'{focal_pop}.vs.Zam_FST',
                    f'{focal_pop}.vs.OOA-OW_FST']


    #removing sites not mapped to a gene
    dfgenes = df[df['FBgn'].notna()]


    #clipping Fst values to zero
    dfgenes.loc[:, comparison_cols] = dfgenes.loc[:, comparison_cols].clip(lower=0)


    #expanding sites with mutliple genes
    dfgenes.loc[:, 'FBgn'] = dfgenes['FBgn'].str.split('; ')
    dfgenes.loc[:, 'gene_symbol'] = dfgenes['gene_symbol'].str.split('; ')
    dfgenes_expand = dfgenes.explode(['FBgn', 'gene_symbol'])

    #average across comparisons
    dfgenes_expand[focal_pop + ".vs.COS"] = dfgenes_expand[comparison_cols].mean(axis=1)



    #averaging Fst per gene
    dfgenes_expand = dfgenes_expand.drop(columns=['CHROM', 'POS'])
    if sumstat == "mean":
        fst_per_gene = dfgenes_expand.groupby(['FBgn', 'gene_symbol']).mean().reset_index()
    elif sumstat == "max":
        fst_per_gene = dfgenes_expand.groupby(['FBgn', 'gene_symbol']).max().reset_index()
    else:
        print("\n\nincorrect sumstat argument\n\n")


    #outputting dfs
    Zim_Cos_gene_Fst = fst_per_gene.drop(columns=comparison_cols)
    Zim_Cos_gene_Fst = Zim_Cos_gene_Fst.sort_values(by=[focal_pop + ".vs.COS"], ascending=False)


    return Zim_Cos_gene_Fst






fstZim_mean = FST_per_gene_across_comparisons("genes_FST/Zim_fst_genes.csv", "Zim", "mean")
fstZS_mean = FST_per_gene_across_comparisons("genes_FST/ZS_fst_genes.csv", "ZS", "mean")
fstZH_mean = FST_per_gene_across_comparisons("genes_FST/ZH_fst_genes.csv", "ZH", "mean")
fstZC_mean = FST_per_gene_across_comparisons("genes_FST/ZC_fst_genes.csv", "ZC", "mean")
fstZR_mean = FST_per_gene_across_comparisons("genes_FST/ZR_fst_genes.csv", "ZR", "mean")


fstZim_max = FST_per_gene_across_comparisons("genes_FST/Zim_fst_genes.csv", "Zim", "max")
fstZS_max = FST_per_gene_across_comparisons("genes_FST/ZS_fst_genes.csv", "ZS", "max")
fstZH_max = FST_per_gene_across_comparisons("genes_FST/ZH_fst_genes.csv", "ZH", "max")
fstZC_max = FST_per_gene_across_comparisons("genes_FST/ZC_fst_genes.csv", "ZC", "max")
fstZR_max = FST_per_gene_across_comparisons("genes_FST/ZR_fst_genes.csv", "ZR", "max")


fstZim_mean.to_csv("Zim.vs.COS_avgfst_ChrX_genes_ranked.csv", index=False)
fstZS_mean.to_csv("ZS.vs.COS_avgfst_ChrX_genes_ranked.csv", index=False)
fstZH_mean.to_csv("ZH.vs.COS_avgfst_ChrX_genes_ranked.csv", index=False)
fstZC_mean.to_csv("ZC.vs.COS_avgfst_ChrX_genes_ranked.csv", index=False)
fstZR_mean.to_csv("ZR.vs.COS_avgfst_ChrX_genes_ranked.csv", index=False)


fstZim_max.to_csv("Zim.vs.COS_maxfst_ChrX_genes_ranked.csv", index=False)
fstZS_max.to_csv("ZS.vs.COS_maxfst_ChrX_genes_ranked.csv", index=False)
fstZH_max.to_csv("ZH.vs.COS_maxfst_ChrX_genes_ranked.csv", index=False)
fstZC_max.to_csv("ZC.vs.COS_maxfst_ChrX_genes_ranked.csv", index=False)
fstZR_max.to_csv("ZR.vs.COS_maxfst_ChrX_genes_ranked.csv", index=False)
