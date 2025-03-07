import pandas as pd



def filter4NS_and_intron_variants(snpEff_genes_file):
    """
    Function takes in the genes tsv file from snpEff and saves a csv for NS variant genes and intro variant genes
    """

    df = pd.read_csv(snpEff_genes_file, sep="\t")

    df_NS = df[(df['variants_effect_missense_variant'] > 0) | (df['variants_effect_stop_gained'] > 0)]

    df_intron = df[df['variants_effect_intron_variant'] > 0]

    df_NS.to_csv(snpEff_genes_file.split("_")[0] + "_NS_variant_genes.csv", index=False)
    df_intron.to_csv(snpEff_genes_file.split("_")[0] + "_intron_variant_genes.csv", index=False)



filter4NS_and_intron_variants("Zim_snpEff_genes.txt")
filter4NS_and_intron_variants("ZS_snpEff_genes.txt")
filter4NS_and_intron_variants("ZH_snpEff_genes.txt")
filter4NS_and_intron_variants("ZC_snpEff_genes.txt")
filter4NS_and_intron_variants("ZR_snpEff_genes.txt")
filter4NS_and_intron_variants("all_snpEff_genes.txt")