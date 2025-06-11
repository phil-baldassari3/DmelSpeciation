import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from statistics import mean


def annotate_genes_with_GO(data, GO_dictionary):
    """
    Function takes in a dataframe of genomic statistics data mapped to genes and annotates sites/genes with
    the Gene Ontology they belong to. The function adds a column for each GO term with values of True or
    False to denote whether the gene belongs to the GO category.

    Arguments:
    data (DataFrame): genomic statistics data
    GO_dictionary (dict): dictionary of GO term keys and gene list values

    Returns Dataframe
    """

    #copying df
    df = data.copy()

    #renaming columns
    df.columns = [col.split("_")[0] for col in df.columns]

    #extracting genes
    genes = df["FBgn"].to_list()

    #annotating genes
    annotations = {}
    for term, term_genes in GO_dictionary.items():
        annotations[term] = []

        for mappedgenes in genes:
            for gene in mappedgenes.split("; "):
                annot = gene in term_genes
                if annot:
                    break
            annotations[term].append(annot)

    #creating annotated df
    annotations_df = pd.DataFrame(annotations)
    df = pd.concat([df, annotations_df], axis=1)

    return df



def avg_stat_perGO_perPop(annotdf, GOtermlist, poplist, sumstat=mean, subset_top=None, subset_bottom=None):
    """
    Function takes in the annotated df and using a list of GO terms and a list of population/comparision
    names (for ordering purposes) outputs a table that can be used to generate a heatmap of average stat
    values per GO per pop.

    Arguments:
    annotdf (DataFrame): annotated df from annotate_genes_with_GO
    GOtermlist (list): list of GO terms used for filering
    poplist (list): list of population or p opulation comparison names used to keep the output df ordered

    Optional Arguments:
    sumstat (function): either mean (from statistics module) or max
    subset_top (int): subset the data for the top n values before comuting sumstat. Cannot be used with subset_bottom.
    subset_bottom (int): subset the data for the bottom n values before comuting sumstat. Cannot be used with subset_top.

    Returns DataFrame
    """

    #filling in table
    table_dict = {}

    for term in GOtermlist:
        col_list = []

        #filtering
        filtered_df = annotdf[annotdf[term] == True]

        for pop in poplist:
            values = filtered_df[pop].to_list()
            
            if subset_top is not None:
                values.sort(reverse=True)
                avg = sumstat(values[:subset_top])
                col_list.append(avg)

            elif subset_bottom is not None:
                values.sort()
                avg = sumstat(values[:subset_bottom])
                col_list.append(avg)

            else:
                avg = sumstat(values)
                col_list.append(avg)


        table_dict[term] = col_list

    #making df
    avgs_df = pd.DataFrame(table_dict)
    avgs_df.insert(0, "", poplist)


    return avgs_df


def plot_heatmap(table, title, stat):
    """
    function plots a heatmap from the table generated from avg_stat_perGO_perPop.
    """

    #setting row labels to index
    table = table.set_index(table.columns[0])
    table.columns = table.columns.str.replace(" system development", "", regex=False)
    table.columns = table.columns.str.replace(" system process", "", regex=False)

    #setting color map
    cmap = plt.cm.plasma

    #plot heat
    plt.figure(figsize=(10, 10))
    heat = plt.imshow(table.values, cmap=cmap)

    #adding scale bar
    cbar = plt.colorbar(heat, shrink=0.4, pad=0.02)
    cbar.set_label(f"Mean {stat}", rotation=270, labelpad=15)

    #adding gridlines
    num_rows, num_cols = table.shape
    for r in range(1, num_rows):
        plt.axhline(r - 0.5, color="white", linewidth=1)
    for c in range(1, num_cols):
        plt.axvline(c - 0.5, color="white", linewidth=1)

    #adding labels
    plt.xticks(range(table.shape[1]), table.columns, rotation=45, ha="right")
    plt.yticks(range(table.shape[0]), table.index)

    #formatting plots
    plt.grid(False)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.title(title)
    plt.subplots_adjust(left=0.2, bottom=0.2)
    #plt.show()

    newtitle = title.replace("\n", "_")
    plt.savefig(f"{newtitle}.png", dpi=300, bbox_inches="tight")
    plt.close()




def main_func(title, genomic_stat, csvfile, GO_directory, popls, summary_statistic=mean, top_subset=None, bottom_subset=None):

    #opening data table
    df = pd.read_csv(csvfile, keep_default_na=False)

    #creating dictioanry of GO terms and genes
    GOdictionary = {}
    for filename in os.listdir(GO_directory):
            if filename.endswith(".txt"):
                file = os.path.join(GO_directory, filename)
                GOterm = " ".join(filename.split("_")).replace("genes.txt", "")
                with open(file, "r") as f:
                    GOgenelist = [line.strip() for line in f.readlines()]
                    GOdictionary[GOterm] = GOgenelist

    #annotating df, extracting data, plotting heatmap
    annotated_df = annotate_genes_with_GO(df, GOdictionary)
    heat_table = avg_stat_perGO_perPop(annotated_df, list(GOdictionary.keys()), popls, sumstat=summary_statistic, subset_top=top_subset, subset_bottom=bottom_subset)
    plot_heatmap(heat_table, title, genomic_stat)




population_list = ["ZS", "ZH", "RAL", "FR", "ZI"]
pop_comparision_list = ["ZS.vs.RAL", "ZS.vs.FR", "ZS.vs.ZI", "ZH.vs.RAL", "ZH.vs.FR", "ZH.vs.ZI", "RAL.vs.FR", "FR.vs.ZI", "ZI.vs.RAL"]


#base paths
fstpath = "FST/FST_data"
pipath = "Pi/Pi_data"
thetapath = "Theta/Theta_data"
tjDpath = "TajimasD/TajimasD_data"



##MEAN

#FST

#System Development
""" main_func(
    "avgFST_snpavg_System_Development", 
    f"{fstpath}/selected_pops_per_gene_fst.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development",
    pop_comparision_list
) """

main_func(
    "FST (Gene-level)\nSystem Development", "FST",
    f"{fstpath}/selected_pops_per_gene_fst_NormByGeneLength.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development",
    pop_comparision_list
)

main_func(
    "FST (10kbp Windows)\nSystem Development", "FST",
    f"{fstpath}/windowed_mean_10kbp_selected_pops_fst_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development",
    pop_comparision_list
)

#Behavior
""" main_func(
    "avgFST_snpavg_Behavior", 
    f"{fstpath}/selected_pops_per_gene_fst.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior",
    pop_comparision_list
) """

main_func(
    "FST (Gene-level)\nBehavior", "FST",
    f"{fstpath}/selected_pops_per_gene_fst_NormByGeneLength.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior",
    pop_comparision_list
)

main_func(
    "FST (10kbp Windows)\nBehavior", "FST",
    f"{fstpath}/windowed_mean_10kbp_selected_pops_fst_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior",
    pop_comparision_list
)


""" #Mating Behavior
main_func(
    "avgFST_snpavg_Mating_Behavior", 
    f"{fstpath}/selected_pops_per_gene_fst.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/mating_behavior",
    pop_comparision_list
)

main_func(
    "avgFST_pergene_Mating_Behavior", 
    f"{fstpath}/selected_pops_per_gene_fst_NormByGeneLength.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/mating_behavior",
    pop_comparision_list
)

main_func(
    "avgFST_windowed_Mating_Behavior", 
    f"{fstpath}/windowed_mean_10kbp_selected_pops_fst_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/mating_behavior",
    pop_comparision_list
) """


#System Process
""" main_func(
    "avgFST_snpavg_System_Process", 
    f"{fstpath}/selected_pops_per_gene_fst.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_process",
    pop_comparision_list
) """

main_func(
    "FST (Gene-level)\nSystem Process", "FST",
    f"{fstpath}/selected_pops_per_gene_fst_NormByGeneLength.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_process",
    pop_comparision_list
)

main_func(
    "FST (10kbp Windows)\nSystem Process", "FST",
    f"{fstpath}/windowed_mean_10kbp_selected_pops_fst_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_process",
    pop_comparision_list
)






#PI

#System Development
main_func(
    "Pi (Gene-level)\nSystem Development", "Pi",
    f"{pipath}/selected_pops_per_gene_pi.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development",
    population_list
)

main_func(
    "Pi (10kbp Windows)\nSystem Development", "Pi",
    f"{pipath}/selected_pops_windowed_10kbp_pi_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development",
    population_list
)

#Behavior
main_func(
    "Pi (Gene-level)\nBehavior", "Pi",
    f"{pipath}/selected_pops_per_gene_pi.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior",
    population_list
)

main_func(
    "Pi (10kbp Windows)\nBehavior", "Pi",
    f"{pipath}/selected_pops_windowed_10kbp_pi_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior",
    population_list
)


""" #Mating Behavior
main_func(
    "avgPI_pergene_Mating_Behavior", 
    f"{pipath}/selected_pops_per_gene_pi.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/mating_behavior",
    population_list
)

main_func(
    "avgPI_windowed_Mating_Behavior", 
    f"{pipath}/selected_pops_windowed_10kbp_pi_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/mating_behavior",
    population_list
) """


#System Process
main_func(
    "Pi (Gene-level)\nSystem Process", "Pi",
    f"{pipath}/selected_pops_per_gene_pi.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_process",
    population_list
)

main_func(
    "Pi (10kbp Windows)\nSystem Process", "Pi",
    f"{pipath}/selected_pops_windowed_10kbp_pi_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_process",
    population_list
)





#Theta

#System Development
main_func(
    "Watterson's Theta (10kbp Windows)\nSystem Development", "Theta",
    f"{thetapath}/selected_pops_windowed_10kbp_theta_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development",
    population_list
)

#Behavior
main_func(
    "Watterson's Theta (10kbp Windows)\nBehavior", "Theta",
    f"{thetapath}/selected_pops_windowed_10kbp_theta_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior",
    population_list
)


""" #Mating Behavior
main_func(
    "avgTheta_Mating_Behavior", 
    f"{thetapath}/selected_pops_windowed_10kbp_theta_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/mating_behavior",
    population_list
) """


#System Process
main_func(
    "Watterson's Theta (10kbp Windows)\nSystem Process", "Theta",
    f"{thetapath}/selected_pops_windowed_10kbp_theta_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_process",
    population_list
)





#TjD

#System Development
main_func(
    "Tajima's D (10kbp Windows)\nSystem_Development", "Tajima's D",
    f"{tjDpath}/selected_pops_windowed_10kbp_TajimasD_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development",
    population_list
)

#Behavior
main_func(
    "Tajima's D (10kbp Windows)\nBehavior", "Tajima's D",
    f"{tjDpath}/selected_pops_windowed_10kbp_TajimasD_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior",
    population_list
)


""" #Mating Behavior
main_func(
    "avgTajimasD_Mating_Behavior", 
    f"{tjDpath}/selected_pops_windowed_10kbp_TajimasD_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/mating_behavior",
    population_list
) """


#System Process
main_func(
    "Tajima's D (10kbp Windows)\nSystem Process", "Tajima's D",
    f"{tjDpath}/selected_pops_windowed_10kbp_TajimasD_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_process",
    population_list
)







##MEAN of subset

#FST

#System Development
""" main_func(
    "top10avgFST_snpavg_System_Development", 
    f"{fstpath}/selected_pops_per_gene_fst.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development",
    pop_comparision_list,
    top_subset=10
) """

main_func(
    "Top 10 FST (Gene-level)\nSystem Development", "of top 10 FST",
    f"{fstpath}/selected_pops_per_gene_fst_NormByGeneLength.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development",
    pop_comparision_list,
    top_subset=10
)

main_func(
    "Top 10 FST (10kbp Windows)\nSystem Development", "of top 10 FST",
    f"{fstpath}/windowed_mean_10kbp_selected_pops_fst_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development",
    pop_comparision_list,
    top_subset=10
)


#Behavior
""" main_func(
    "top10avgFST_snpavg_Behavior", 
    f"{fstpath}/selected_pops_per_gene_fst.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior",
    pop_comparision_list, 
    top_subset=10
) """

main_func(
    "Top 10 FST (Gene-level)\nBehavior", "of top 10 FST",
    f"{fstpath}/selected_pops_per_gene_fst_NormByGeneLength.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior",
    pop_comparision_list, 
    top_subset=10
)

main_func(
    "Top 10 FST (10kbp Windows)\nBehavior", "of top 10 FST",
    f"{fstpath}/windowed_mean_10kbp_selected_pops_fst_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior",
    pop_comparision_list, 
    top_subset=10
)


""" #Mating Behavior
main_func(
    "top10avgFST_snpavg_Mating_Behavior", 
    f"{fstpath}/selected_pops_per_gene_fst.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/mating_behavior",
    pop_comparision_list, 
    top_subset=10
)

main_func(
    "top10avgFST_pergene_Mating_Behavior", 
    f"{fstpath}/selected_pops_per_gene_fst_NormByGeneLength.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/mating_behavior",
    pop_comparision_list, 
    top_subset=10
)

main_func(
    "top10avgFST_windowed_Mating_Behavior", 
    f"{fstpath}/windowed_mean_10kbp_selected_pops_fst_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/mating_behavior",
    pop_comparision_list, 
    top_subset=10
) """


#System Process
""" main_func(
    "top10avgFST_snpavg_System_Process", 
    f"{fstpath}/selected_pops_per_gene_fst.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_process",
    pop_comparision_list,
    top_subset=10
) """

main_func(
    "Top 10 FST (Gene-level)\nSystem Process", "of top 10 FST",
    f"{fstpath}/selected_pops_per_gene_fst_NormByGeneLength.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_process",
    pop_comparision_list,
    top_subset=10
)

main_func(
    "Top 10 FST (10kbp Windows)\nSystem Process", "of top 10 FST",
    f"{fstpath}/windowed_mean_10kbp_selected_pops_fst_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_process",
    pop_comparision_list,
    top_subset=10
)






#PI

#System Development
main_func(
    "Bottom 10 Pi (Gene-level)\nSystem Development", "of bottom 10 Pi",
    f"{pipath}/selected_pops_per_gene_pi.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development",
    population_list,
    bottom_subset=10
    
)

main_func(
    "Bottom 10 Pi (10kbp Windows)\nSystem Development", "of bottom 10 Pi",
    f"{pipath}/selected_pops_windowed_10kbp_pi_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development",
    population_list,
    bottom_subset=10
    
)

#Behavior
main_func(
    "Bottom 10 Pi (Gene-level)\nBehavior", "of bottom 10 Pi",
    f"{pipath}/selected_pops_per_gene_pi.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior",
    population_list,
    bottom_subset=10
    
)

main_func(
    "Bottom 10 Pi (10kbp Windows)\nBehavior", "of bottom 10 Pi",
    f"{pipath}/selected_pops_windowed_10kbp_pi_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior",
    population_list,
    bottom_subset=10
    
)


""" #Mating Behavior
main_func(
    "bottom10avgPI_pergene_Mating_Behavior", 
    f"{pipath}/selected_pops_per_gene_pi.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/mating_behavior",
    population_list,
    bottom_subset=10
    
)

main_func(
    "bottom10avgPI_windowed_Mating_Behavior", 
    f"{pipath}/selected_pops_windowed_10kbp_pi_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/mating_behavior",
    population_list,
    bottom_subset=10
    
) """


#System Process
main_func(
    "Bottom 10 Pi (Gene-level)\nSystem Process", "of bottom 10 Pi",
    f"{pipath}/selected_pops_per_gene_pi.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_process",
    population_list,
    bottom_subset=10
    
)

main_func(
    "Bottom 10 Pi (10kbp Windows)\nSystem Process", "of bottom 10 Pi",
    f"{pipath}/selected_pops_windowed_10kbp_pi_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_process",
    population_list,
    bottom_subset=10
    
)





#Theta

main_func(
    "Bottom 10 Watterson's Theta (10kbp Windows)\nSystem Development", "of bottom 10 Theta",
    f"{thetapath}/selected_pops_windowed_10kbp_theta_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development",
    population_list,
    bottom_subset=10
)

#Behavior
main_func(
    "Bottom 10 Watterson's Theta (10kbp Windows)\nBehavior", "of bottom 10 Theta",
    f"{thetapath}/selected_pops_windowed_10kbp_theta_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior",
    population_list,
    bottom_subset=10
)


""" #Mating Behavior
main_func(
    "bottom10avgTheta_Mating_Behavior", 
    f"{thetapath}/selected_pops_windowed_10kbp_theta_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/mating_behavior",
    population_list,
    bottom_subset=10
) """


#System Process
main_func(
    "Bottom 10 Watterson's Theta (10kbp Windows)\nSystem Process", "of bottom 10 Theta",
    f"{thetapath}/selected_pops_windowed_10kbp_theta_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_process",
    population_list,
    bottom_subset=10
)





#TjD

#System Development
main_func(
    "Bottom 10 Tajima's D (10kbp Windows)\nSystem_Development", "of bottom 10 Tajima's D",
    f"{tjDpath}/selected_pops_windowed_10kbp_TajimasD_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_development",
    population_list,
    bottom_subset=10
)

#Behavior
main_func(
    "Bottom 10 Tajima's D (10kbp Windows)\nBehavior", "of bottom 10 Tajima's D",
    f"{tjDpath}/selected_pops_windowed_10kbp_TajimasD_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/behavior",
    population_list,
    bottom_subset=10
)


""" #Mating Behavior
main_func(
    "bottom10avgTajimasD_Mating_Behavior", 
    f"{tjDpath}/selected_pops_windowed_10kbp_TajimasD_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/mating_behavior",
    population_list,
    bottom_subset=10
) """


#System Process
main_func(
    "Bottom 10 Tajima's D (10kbp Windows)\nSystem Process", "of bottom 10 Tajima's D",
    f"{tjDpath}/selected_pops_windowed_10kbp_TajimasD_mapped2genes.csv", 
    "GO_enrichment_and_annotation/lists_of_genes_per_GO_category_from_flymine/system_process",
    population_list,
    bottom_subset=10
)