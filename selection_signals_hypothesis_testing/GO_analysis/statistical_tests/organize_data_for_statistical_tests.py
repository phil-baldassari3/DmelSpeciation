import pandas as pd



def annotate_genes_with_GO(df, GO_dictionary):
    """
    Function takes in a dataframe of genomic statistics data mapped to genes and annotates sites/genes with
    the Gene Ontology they belong to. The function adds a column for each GO term with values of "term" or
    "non-term" to denote whether the gene belongs to the GO category or not. The function also filters the
    dataframe for rows mapped to genes only.

    Arguments:
    data (DataFrame): genomic statistics data
    GO_dictionary (dict): dictionary of GO term keys and gene list values

    Returns Dataframe
    """

    #renaming columns
    df.columns = [col.split("_")[0] for col in df.columns]

    #filtering for rows mapped to genes only
    df = df[df["FBgn"] != "None"]

    #extracting genes
    genes = df["FBgn"].to_list()

    #annotating genes
    annotations = {}
    for term, term_genes in GO_dictionary.items():
        annotations[term] = []

        for mappedgenes in genes:
            for gene in mappedgenes.split("; "):
                if gene in term_genes:
                    annot = " ".join(term.split("_"))
                    break
                else:
                    annot = f"Non-" + " ".join(term.split("_"))
            annotations[term].append(annot)

    #creating annotated df
    df = df.reset_index(drop=True)
    annotations_df = pd.DataFrame(annotations)
    df = pd.concat([df, annotations_df], axis=1)

    #dropping unused columns
    if "CHROM" in df.columns:
        df.drop(columns=["CHROM", "START", "END", "FBgn", "gene"], inplace=True)
    else:
        df.drop(columns=["FBgn", "gene"], inplace=True)

    return df






def main_func(csvfile, GO_files, outfilename):

    #opening data table
    df = pd.read_csv(csvfile, keep_default_na=False)

    #creating dictioanry of GO terms and genes
    GOdictionary = {}
    for filename in GO_files:
        GOterm = filename.split("/")[-1].replace("_genes.txt", "")
        with open(filename, "r") as f:
            GOgenelist = [line.strip() for line in f.readlines()]
            GOdictionary[GOterm] = GOgenelist

    #annotating df, extracting data, plotting heatmap
    annotated_df = annotate_genes_with_GO(df, GOdictionary)

    #saving csv
    annotated_df.to_csv(f"GO_ann_{outfilename}.csv", index=False)







#base paths
fstpath = "FST/FST_data"
pipath = "Pi/Pi_data"
thetapath = "Theta/Theta_data"
tjDpath = "TajimasD/TajimasD_data"


#GO directories
GOs = [
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development/nervous_system_development_genes.txt",
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development/sensory_system_development_genes.txt",
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior/mechanosensory_behavior_genes.txt",
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior/aggressive_behavior_genes.txt",
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior/learning_or_memory_genes.txt",
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_process/nervous_system_process_genes.txt"
]


main_func(f"{fstpath}/selected_pops_per_gene_fst_NormByGeneLength.csv", GOs, "FST_gene")
main_func(f"{fstpath}/windowed_mean_10kbp_selected_pops_fst_mapped2genes.csv", GOs, "FST_win")
main_func(f"{pipath}/selected_pops_per_gene_pi.csv", GOs, "PI_gene")
main_func(f"{pipath}/selected_pops_windowed_10kbp_pi_mapped2genes.csv", GOs, "PI_win")
main_func(f"{thetapath}/selected_pops_windowed_10kbp_theta_mapped2genes.csv", GOs, "THETA")
main_func(f"{tjDpath}/selected_pops_windowed_10kbp_TajimasD_mapped2genes.csv", GOs, "TjD")