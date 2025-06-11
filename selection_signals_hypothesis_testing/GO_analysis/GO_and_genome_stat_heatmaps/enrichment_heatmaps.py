import pandas as pd
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np



def enrichment_heatmap(table, title):
    """
    """

    #setting row labels to index
    table = table.set_index(table.columns[0])
    table.columns = table.columns.str.replace(" system development", "", regex=False)

    #setting color map
    base_cmap = plt.cm.Greens
    new_colors = base_cmap(np.linspace(0, 1, 256))
    new_colors[0] = [1, 1, 1, 1]  # override 0 with white
    cmap = mcolors.ListedColormap(new_colors)
    norm = mcolors.Normalize(vmin=0, vmax=10)

    #plot heat
    plt.figure(figsize=(10, 10))
    plt.imshow(table.values, cmap=cmap, norm=norm)

    #adding labels
    plt.xticks(range(table.shape[1]), table.columns, rotation=45, ha="right")
    plt.yticks(range(table.shape[0]), table.index)

    #adding counts
    count_data = table.values.astype(int)
    for i in range(count_data.shape[0]):
        for j in range(count_data.shape[1]):
            val = count_data[i, j]
            plt.text(j, i, str(val), ha="center", va="center", fontsize=10, color="black")

    #formatting plots
    plt.grid(False)
    plt.gca().set_aspect('equal', adjustable='box')

    plt.title(title)
    plt.subplots_adjust(left=0.2, bottom=0.2)
    #plt.show()

    plt.savefig(f"{title}.png", dpi=300, bbox_inches="tight")





def enrichment_counts_table(name, enrichment_results_directory, GO_child_terms_directory, pop_list):
    """
    Function makes a table showing how many child terms (3-levels) a population/population comparison is enriched for
    a specific GO term (3-levels)

    returns DataFrame

    e.g.
          GO1  GO2
    Pop1  5    0
    Pop2  0    0
    Pop3  2    1
    
    Arguments:
    name (str): label for table used for saving csv and title for plotting
    enrichment_results_directory (str): path to directory containing lists (.txt) of enriched GO terms for each population/population comparision
    GO_child_terms_directory (str): path to directory containing lists (.txt) of child terms under each GO term
    pop_list (list): list of population/population comparison labels for ordering purposes. must match pop names in enrichment_results_directory
    """

    #enrichment dictionary
    enrichments = {}
    for filename in os.listdir(enrichment_results_directory):
        if filename.endswith(".txt"):
            file = os.path.join(enrichment_results_directory, filename)

            pop_label = filename.split("_")[0]
            
            with open(file, "r") as f:
                enrichment_list = [line.strip() for line in f.readlines()]

                enrichments[pop_label] = enrichment_list

    #GO child term dictioanry
    go = {}
    for filename in os.listdir(GO_child_terms_directory):
        if filename.endswith(".txt"):
            file = os.path.join(GO_child_terms_directory, filename)

            GOterm = " ".join(filename.split("_")).replace(".txt", "")
            
            with open(file, "r") as f:
                GOchildren_list = [line.strip() for line in f.readlines()]

                go[GOterm] = GOchildren_list

    #GO terms list (for maintaining table order)
    GOterms_list = list(go.keys())
    
    #filling in table
    table_dict = {}

    for term in GOterms_list:

        col_list = []
        for pop in pop_list:

            counter = 0
            for enriched_term in enrichments[pop]:
                if enriched_term in go[term]:
                    counter += 1
                    #break

            col_list.append(counter)

        table_dict[term] = col_list

    #making table
    enrichment_table = pd.DataFrame(table_dict)
    enrichment_table.insert(0, "", pop_list)

    #enrichment_table.to_csv(f"{name}.csv", index=False)
    #enrichment_table.to_excel(f"{name}.xlsx", index=False)


    #plotting
    enrichment_heatmap(enrichment_table, name)



    return enrichment_table









population_list = ["ZS", "ZH", "RAL", "FR", "ZI"]
pop_comparision_list = ["ZS.vs.RAL", "ZS.vs.FR", "ZS.vs.ZI", "ZH.vs.RAL", "ZH.vs.FR", "ZH.vs.ZI", "RAL.vs.FR", "FR.vs.ZI", "ZI.vs.RAL"]




#top/bottom 1% System Development Enrichment
enrichment_counts_table("top1_geneFST_system_development", "top_bottom_1percent_genes/enrichment_results/fst_enriched_GO/pergene_fst", "GO_child_terms/system_development", pop_comparision_list)

enrichment_counts_table("bottom1_genePI_system_development", "top_bottom_1percent_genes/enrichment_results/pi_enriched_GO/pergene_pi", "GO_child_terms/system_development", population_list)
enrichment_counts_table("bottom1_windowPI_system_development", "top_bottom_1percent_genes/enrichment_results/pi_enriched_GO/windowed_pi", "GO_child_terms/system_development", population_list)

enrichment_counts_table("bottom1_THETA_system_development", "top_bottom_1percent_genes/enrichment_results/theta_enriched_GO", "GO_child_terms/system_development", population_list)

enrichment_counts_table("bottom1_TAJIMA_system_development", "top_bottom_1percent_genes/enrichment_results/tajimasD_enriched_GO", "GO_child_terms/system_development", population_list)

#top/bottom 1% Behavior Enrichment
enrichment_counts_table("top1_geneFST_behavior", "top_bottom_1percent_genes/enrichment_results/fst_enriched_GO/pergene_fst", "GO_child_terms/behavior", pop_comparision_list)

enrichment_counts_table("bottom1_genePI_behavior", "top_bottom_1percent_genes/enrichment_results/pi_enriched_GO/pergene_pi", "GO_child_terms/behavior", population_list)
enrichment_counts_table("bottom1_windowPI_behavior", "top_bottom_1percent_genes/enrichment_results/pi_enriched_GO/windowed_pi", "GO_child_terms/behavior", population_list)

enrichment_counts_table("bottom1_THETA_behavior", "top_bottom_1percent_genes/enrichment_results/theta_enriched_GO", "GO_child_terms/behavior", population_list)

enrichment_counts_table("bottom1_TAJIMA_behavior", "top_bottom_1percent_genes/enrichment_results/tajimasD_enriched_GO", "GO_child_terms/behavior", population_list)

#top/bottom 1% Mating behavior Enrichment
enrichment_counts_table("top1_geneFST_mating_behavior", "top_bottom_1percent_genes/enrichment_results/fst_enriched_GO/pergene_fst", "GO_child_terms/mating_behavior", pop_comparision_list)

enrichment_counts_table("bottom1_genePI_mating_behavior", "top_bottom_1percent_genes/enrichment_results/pi_enriched_GO/pergene_pi", "GO_child_terms/mating_behavior", population_list)
enrichment_counts_table("bottom1_windowPI_mating_behavior", "top_bottom_1percent_genes/enrichment_results/pi_enriched_GO/windowed_pi", "GO_child_terms/mating_behavior", population_list)

enrichment_counts_table("bottom1_THETA_mating_behavior", "top_bottom_1percent_genes/enrichment_results/theta_enriched_GO", "GO_child_terms/mating_behavior", population_list)

enrichment_counts_table("bottom1_TAJIMA_mating_behavior", "top_bottom_1percent_genes/enrichment_results/tajimasD_enriched_GO", "GO_child_terms/mating_behavior", population_list)

#top/bottom 1% System Process
enrichment_counts_table("top1_geneFST_system_process", "top_bottom_1percent_genes/enrichment_results/fst_enriched_GO/pergene_fst", "GO_child_terms/system_process", pop_comparision_list)

enrichment_counts_table("bottom1_genePI_system_process", "top_bottom_1percent_genes/enrichment_results/pi_enriched_GO/pergene_pi", "GO_child_terms/system_process", population_list)
enrichment_counts_table("bottom1_windowPI_system_process", "top_bottom_1percent_genes/enrichment_results/pi_enriched_GO/windowed_pi", "GO_child_terms/system_process", population_list)

enrichment_counts_table("bottom1_THETA_system_process", "top_bottom_1percent_genes/enrichment_results/theta_enriched_GO", "GO_child_terms/system_process", population_list)

enrichment_counts_table("bottom1_TAJIMA_system_process", "top_bottom_1percent_genes/enrichment_results/tajimasD_enriched_GO", "GO_child_terms/system_process", population_list)



#top/bottom 5% System Development Enrichment
enrichment_counts_table("top5_geneFST_system_development", "top_bottom_5percent_genes/enrichment_results/fst_enriched_GO/pergene_fst", "GO_child_terms/system_development", pop_comparision_list)
enrichment_counts_table("top5_windowFST_system_development", "top_bottom_5percent_genes/enrichment_results/fst_enriched_GO/windowed_fst", "GO_child_terms/system_development", pop_comparision_list)

enrichment_counts_table("bottom5_genePI_system_development", "top_bottom_5percent_genes/enrichment_results/pi_enriched_GO/pergene_pi", "GO_child_terms/system_development", population_list)
enrichment_counts_table("bottom5_windowPI_system_development", "top_bottom_5percent_genes/enrichment_results/pi_enriched_GO/windowed_pi", "GO_child_terms/system_development", population_list)

enrichment_counts_table("bottom5_THETA_system_development", "top_bottom_5percent_genes/enrichment_results/theta_enriched_GO", "GO_child_terms/system_development", population_list)

enrichment_counts_table("bottom5_TAJIMA_system_development", "top_bottom_5percent_genes/enrichment_results/tajimasD_enriched_GO", "GO_child_terms/system_development", population_list)

#top/bottom 5% Behavior Enrichment
enrichment_counts_table("top5_geneFST_behavior", "top_bottom_5percent_genes/enrichment_results/fst_enriched_GO/pergene_fst", "GO_child_terms/behavior", pop_comparision_list)
enrichment_counts_table("top5_windowFST_behavior", "top_bottom_5percent_genes/enrichment_results/fst_enriched_GO/windowed_fst", "GO_child_terms/behavior", pop_comparision_list)

enrichment_counts_table("bottom5_genePI_behavior", "top_bottom_5percent_genes/enrichment_results/pi_enriched_GO/pergene_pi", "GO_child_terms/behavior", population_list)
enrichment_counts_table("bottom5_windowPI_behavior", "top_bottom_5percent_genes/enrichment_results/pi_enriched_GO/windowed_pi", "GO_child_terms/behavior", population_list)

enrichment_counts_table("bottom5_THETA_behavior", "top_bottom_5percent_genes/enrichment_results/theta_enriched_GO", "GO_child_terms/behavior", population_list)

enrichment_counts_table("bottom5_TAJIMA_behavior", "top_bottom_5percent_genes/enrichment_results/tajimasD_enriched_GO", "GO_child_terms/behavior", population_list)

#top/bottom 5% Mating behavior Enrichment
enrichment_counts_table("top5_geneFST_mating_behavior", "top_bottom_5percent_genes/enrichment_results/fst_enriched_GO/pergene_fst", "GO_child_terms/mating_behavior", pop_comparision_list)
enrichment_counts_table("top5_windowFST_mating_behavior", "top_bottom_5percent_genes/enrichment_results/fst_enriched_GO/windowed_fst", "GO_child_terms/mating_behavior", pop_comparision_list)

enrichment_counts_table("bottom5_genePI_mating_behavior", "top_bottom_5percent_genes/enrichment_results/pi_enriched_GO/pergene_pi", "GO_child_terms/mating_behavior", population_list)
enrichment_counts_table("bottom5_windowPI_mating_behavior", "top_bottom_5percent_genes/enrichment_results/pi_enriched_GO/windowed_pi", "GO_child_terms/mating_behavior", population_list)

enrichment_counts_table("bottom5_THETA_mating_behavior", "top_bottom_5percent_genes/enrichment_results/theta_enriched_GO", "GO_child_terms/mating_behavior", population_list)

enrichment_counts_table("bottom5_TAJIMA_mating_behavior", "top_bottom_5percent_genes/enrichment_results/tajimasD_enriched_GO", "GO_child_terms/mating_behavior", population_list)

#top/bottom 5% System Process
enrichment_counts_table("top5_geneFST_system_process", "top_bottom_5percent_genes/enrichment_results/fst_enriched_GO/pergene_fst", "GO_child_terms/system_process", pop_comparision_list)
enrichment_counts_table("top5_windowFST_system_process", "top_bottom_5percent_genes/enrichment_results/fst_enriched_GO/windowed_fst", "GO_child_terms/system_process", pop_comparision_list)

enrichment_counts_table("bottom5_genePI_system_process", "top_bottom_5percent_genes/enrichment_results/pi_enriched_GO/pergene_pi", "GO_child_terms/system_process", population_list)
enrichment_counts_table("bottom5_windowPI_system_process", "top_bottom_5percent_genes/enrichment_results/pi_enriched_GO/windowed_pi", "GO_child_terms/system_process", population_list)

enrichment_counts_table("bottom5_THETA_system_process", "top_bottom_5percent_genes/enrichment_results/theta_enriched_GO", "GO_child_terms/system_process", population_list)

enrichment_counts_table("bottom5_TAJIMA_system_process", "top_bottom_5percent_genes/enrichment_results/tajimasD_enriched_GO", "GO_child_terms/system_process", population_list)





#top/bottom 10% System Development Enrichment
enrichment_counts_table("top10_geneFST_system_development", "top_bottom_10percent_genes/enrichment_results/fst_enriched_GO/pergene_fst", "GO_child_terms/system_development", pop_comparision_list)
enrichment_counts_table("top10_windowFST_system_development", "top_bottom_10percent_genes/enrichment_results/fst_enriched_GO/windowed_fst", "GO_child_terms/system_development", pop_comparision_list)

enrichment_counts_table("bottom10_genePI_system_development", "top_bottom_10percent_genes/enrichment_results/pi_enriched_GO/pergene_pi", "GO_child_terms/system_development", population_list)
enrichment_counts_table("bottom10_windowPI_system_development", "top_bottom_10percent_genes/enrichment_results/pi_enriched_GO/windowed_pi", "GO_child_terms/system_development", population_list)

enrichment_counts_table("bottom10_THETA_system_development", "top_bottom_10percent_genes/enrichment_results/theta_enriched_GO", "GO_child_terms/system_development", population_list)

enrichment_counts_table("bottom10_TAJIMA_system_development", "top_bottom_10percent_genes/enrichment_results/tajimasD_enriched_GO", "GO_child_terms/system_development", population_list)

#top/bottom 10% Behavior Enrichment
enrichment_counts_table("top10_geneFST_behavior", "top_bottom_10percent_genes/enrichment_results/fst_enriched_GO/pergene_fst", "GO_child_terms/behavior", pop_comparision_list)
enrichment_counts_table("top10_windowFST_behavior", "top_bottom_10percent_genes/enrichment_results/fst_enriched_GO/windowed_fst", "GO_child_terms/behavior", pop_comparision_list)

enrichment_counts_table("bottom10_genePI_behavior", "top_bottom_10percent_genes/enrichment_results/pi_enriched_GO/pergene_pi", "GO_child_terms/behavior", population_list)
enrichment_counts_table("bottom10_windowPI_behavior", "top_bottom_10percent_genes/enrichment_results/pi_enriched_GO/windowed_pi", "GO_child_terms/behavior", population_list)

enrichment_counts_table("bottom10_THETA_behavior", "top_bottom_10percent_genes/enrichment_results/theta_enriched_GO", "GO_child_terms/behavior", population_list)

enrichment_counts_table("bottom10_TAJIMA_behavior", "top_bottom_10percent_genes/enrichment_results/tajimasD_enriched_GO", "GO_child_terms/behavior", population_list)

#top/bottom 10% Mating behavior Enrichment
enrichment_counts_table("top10_geneFST_mating_behavior", "top_bottom_10percent_genes/enrichment_results/fst_enriched_GO/pergene_fst", "GO_child_terms/mating_behavior", pop_comparision_list)
enrichment_counts_table("top10_windowFST_mating_behavior", "top_bottom_10percent_genes/enrichment_results/fst_enriched_GO/windowed_fst", "GO_child_terms/mating_behavior", pop_comparision_list)

enrichment_counts_table("bottom10_genePI_mating_behavior", "top_bottom_10percent_genes/enrichment_results/pi_enriched_GO/pergene_pi", "GO_child_terms/mating_behavior", population_list)
enrichment_counts_table("bottom10_windowPI_mating_behavior", "top_bottom_10percent_genes/enrichment_results/pi_enriched_GO/windowed_pi", "GO_child_terms/mating_behavior", population_list)

enrichment_counts_table("bottom10_THETA_mating_behavior", "top_bottom_10percent_genes/enrichment_results/theta_enriched_GO", "GO_child_terms/mating_behavior", population_list)

enrichment_counts_table("bottom10_TAJIMA_mating_behavior", "top_bottom_10percent_genes/enrichment_results/tajimasD_enriched_GO", "GO_child_terms/mating_behavior", population_list)

#top/bottom 10% System Process
enrichment_counts_table("top10_geneFST_system_process", "top_bottom_10percent_genes/enrichment_results/fst_enriched_GO/pergene_fst", "GO_child_terms/system_process", pop_comparision_list)
enrichment_counts_table("top10_windowFST_system_process", "top_bottom_10percent_genes/enrichment_results/fst_enriched_GO/windowed_fst", "GO_child_terms/system_process", pop_comparision_list)

enrichment_counts_table("bottom10_genePI_system_process", "top_bottom_10percent_genes/enrichment_results/pi_enriched_GO/pergene_pi", "GO_child_terms/system_process", population_list)
enrichment_counts_table("bottom10_windowPI_system_process", "top_bottom_10percent_genes/enrichment_results/pi_enriched_GO/windowed_pi", "GO_child_terms/system_process", population_list)

enrichment_counts_table("bottom10_THETA_system_process", "top_bottom_10percent_genes/enrichment_results/theta_enriched_GO", "GO_child_terms/system_process", population_list)

enrichment_counts_table("bottom10_TAJIMA_system_process", "top_bottom_10percent_genes/enrichment_results/tajimasD_enriched_GO", "GO_child_terms/system_process", population_list)