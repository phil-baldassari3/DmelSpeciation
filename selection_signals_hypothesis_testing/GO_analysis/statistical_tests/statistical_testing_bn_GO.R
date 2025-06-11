library(tidyverse)
library(ggplot2)
library(car)
library(stringr)
library(dplyr)
library(tidyr)

#setting wd
setwd("~/Desktop/dros_speciation/selection_signals_hypothesis_testing/statistical_testing_between_GO_categories")

#Loading Data
fst_gene <- read.csv("GO_ann_FST_gene.csv")
fst_win <- read.csv("GO_ann_FST_win.csv")
pi_gene <- read.csv("GO_ann_PI_gene.csv")
pi_win <- read.csv("GO_ann_PI_win.csv")
theta <- read.csv("GO_ann_THETA.csv")
tajimaD <- read.csv("GO_ann_TjD.csv")


#function parses data, grabs the 2 groups, and performs the statistical test
run_wilcoxon_test <- function(data, GOterm, popname, direction = "two.sided"){
  
  GOterm_value <- str_replace_all(GOterm, "_", " ")
  
  test_group <- data[[popname]][data[[GOterm]] == GOterm_value]
  control_group <- data[[popname]][data[[GOterm]] != GOterm_value]
  
  #print(wilcox.test(test_group, control_group, alternative = "two.sided"))
  print(wilcox.test(test_group, control_group, alternative = direction))
  
}




########running tests########################################################################

#FST_gene#################
#nerv sys dev
run_wilcoxon_test(fst_gene, "nervous_system_development", "ZS.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_development", "ZS.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_development", "ZS.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_development", "ZH.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_development", "ZH.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_development", "ZH.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_development", "RAL.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_development", "FR.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_development", "ZI.vs.RAL", direction = "greater")

#sensory sys dev
run_wilcoxon_test(fst_gene, "sensory_system_development", "ZS.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_gene, "sensory_system_development", "ZS.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "sensory_system_development", "ZS.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "sensory_system_development", "ZH.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_gene, "sensory_system_development", "ZH.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "sensory_system_development", "ZH.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "sensory_system_development", "RAL.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "sensory_system_development", "FR.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "sensory_system_development", "ZI.vs.RAL", direction = "greater")

#mechanosensory behavior
run_wilcoxon_test(fst_gene, "mechanosensory_behavior", "ZS.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_gene, "mechanosensory_behavior", "ZS.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "mechanosensory_behavior", "ZS.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "mechanosensory_behavior", "ZH.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_gene, "mechanosensory_behavior", "ZH.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "mechanosensory_behavior", "ZH.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "mechanosensory_behavior", "RAL.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "mechanosensory_behavior", "FR.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "mechanosensory_behavior", "ZI.vs.RAL", direction = "greater")

#aggressive behavior
run_wilcoxon_test(fst_gene, "aggressive_behavior", "ZS.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_gene, "aggressive_behavior", "ZS.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "aggressive_behavior", "ZS.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "aggressive_behavior", "ZH.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_gene, "aggressive_behavior", "ZH.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "aggressive_behavior", "ZH.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "aggressive_behavior", "RAL.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "aggressive_behavior", "FR.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "aggressive_behavior", "ZI.vs.RAL", direction = "greater")

#learning or memory
run_wilcoxon_test(fst_gene, "learning_or_memory", "ZS.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_gene, "learning_or_memory", "ZS.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "learning_or_memory", "ZS.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "learning_or_memory", "ZH.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_gene, "learning_or_memory", "ZH.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "learning_or_memory", "ZH.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "learning_or_memory", "RAL.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "learning_or_memory", "FR.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "learning_or_memory", "ZI.vs.RAL", direction = "greater")

#nerv sys proc
run_wilcoxon_test(fst_gene, "nervous_system_process", "ZS.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_process", "ZS.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_process", "ZS.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_process", "ZH.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_process", "ZH.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_process", "ZH.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_process", "RAL.vs.FR", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_process", "FR.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_gene, "nervous_system_process", "ZI.vs.RAL", direction = "greater")


###FST_win#################
#nerv sys dev
run_wilcoxon_test(fst_win, "nervous_system_development", "ZS.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_development", "ZS.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_development", "ZS.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_development", "ZH.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_development", "ZH.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_development", "ZH.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_development", "RAL.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_development", "FR.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_development", "ZI.vs.RAL", direction = "greater")

#sensory sys dev
run_wilcoxon_test(fst_win, "sensory_system_development", "ZS.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_win, "sensory_system_development", "ZS.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "sensory_system_development", "ZS.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "sensory_system_development", "ZH.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_win, "sensory_system_development", "ZH.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "sensory_system_development", "ZH.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "sensory_system_development", "RAL.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "sensory_system_development", "FR.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "sensory_system_development", "ZI.vs.RAL", direction = "greater")

#mechanosensory behavior
run_wilcoxon_test(fst_win, "mechanosensory_behavior", "ZS.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_win, "mechanosensory_behavior", "ZS.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "mechanosensory_behavior", "ZS.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "mechanosensory_behavior", "ZH.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_win, "mechanosensory_behavior", "ZH.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "mechanosensory_behavior", "ZH.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "mechanosensory_behavior", "RAL.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "mechanosensory_behavior", "FR.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "mechanosensory_behavior", "ZI.vs.RAL", direction = "greater")

#aggressive behavior
run_wilcoxon_test(fst_win, "aggressive_behavior", "ZS.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_win, "aggressive_behavior", "ZS.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "aggressive_behavior", "ZS.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "aggressive_behavior", "ZH.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_win, "aggressive_behavior", "ZH.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "aggressive_behavior", "ZH.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "aggressive_behavior", "RAL.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "aggressive_behavior", "FR.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "aggressive_behavior", "ZI.vs.RAL", direction = "greater")

#learning or memory
run_wilcoxon_test(fst_win, "learning_or_memory", "ZS.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_win, "learning_or_memory", "ZS.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "learning_or_memory", "ZS.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "learning_or_memory", "ZH.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_win, "learning_or_memory", "ZH.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "learning_or_memory", "ZH.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "learning_or_memory", "RAL.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "learning_or_memory", "FR.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "learning_or_memory", "ZI.vs.RAL", direction = "greater")

#nerv sys proc
run_wilcoxon_test(fst_win, "nervous_system_process", "ZS.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_process", "ZS.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_process", "ZS.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_process", "ZH.vs.RAL", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_process", "ZH.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_process", "ZH.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_process", "RAL.vs.FR", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_process", "FR.vs.ZI", direction = "greater")
run_wilcoxon_test(fst_win, "nervous_system_process", "ZI.vs.RAL", direction = "greater")


###PI_gene#################
#nerv sys dev
run_wilcoxon_test(pi_gene, "nervous_system_development", "ZS", direction = "less")
run_wilcoxon_test(pi_gene, "nervous_system_development", "ZH", direction = "less")
run_wilcoxon_test(pi_gene, "nervous_system_development", "RAL", direction = "less")
run_wilcoxon_test(pi_gene, "nervous_system_development", "FR", direction = "less")
run_wilcoxon_test(pi_gene, "nervous_system_development", "ZI", direction = "less")

#sensory sys dev
run_wilcoxon_test(pi_gene, "sensory_system_development", "ZS", direction = "less")
run_wilcoxon_test(pi_gene, "sensory_system_development", "ZH", direction = "less")
run_wilcoxon_test(pi_gene, "sensory_system_development", "RAL", direction = "less")
run_wilcoxon_test(pi_gene, "sensory_system_development", "FR", direction = "less")
run_wilcoxon_test(pi_gene, "sensory_system_development", "ZI", direction = "less")

#mechanosensory behavior
run_wilcoxon_test(pi_gene, "mechanosensory_behavior", "ZS", direction = "less")
run_wilcoxon_test(pi_gene, "mechanosensory_behavior", "ZH", direction = "less")
run_wilcoxon_test(pi_gene, "mechanosensory_behavior", "RAL", direction = "less")
run_wilcoxon_test(pi_gene, "mechanosensory_behavior", "FR", direction = "less")
run_wilcoxon_test(pi_gene, "mechanosensory_behavior", "ZI", direction = "less")

#aggressive behavior
run_wilcoxon_test(pi_gene, "aggressive_behavior", "ZS", direction = "less")
run_wilcoxon_test(pi_gene, "aggressive_behavior", "ZH", direction = "less")
run_wilcoxon_test(pi_gene, "aggressive_behavior", "RAL", direction = "less")
run_wilcoxon_test(pi_gene, "aggressive_behavior", "FR", direction = "less")
run_wilcoxon_test(pi_gene, "aggressive_behavior", "ZI", direction = "less")

#learning or memory
run_wilcoxon_test(pi_gene, "learning_or_memory", "ZS", direction = "less")
run_wilcoxon_test(pi_gene, "learning_or_memory", "ZH", direction = "less")
run_wilcoxon_test(pi_gene, "learning_or_memory", "RAL", direction = "less")
run_wilcoxon_test(pi_gene, "learning_or_memory", "FR", direction = "less")
run_wilcoxon_test(pi_gene, "learning_or_memory", "ZI", direction = "less")

#nerv sys proc
run_wilcoxon_test(pi_gene, "nervous_system_process", "ZS", direction = "less")
run_wilcoxon_test(pi_gene, "nervous_system_process", "ZH", direction = "less")
run_wilcoxon_test(pi_gene, "nervous_system_process", "RAL", direction = "less")
run_wilcoxon_test(pi_gene, "nervous_system_process", "FR", direction = "less")
run_wilcoxon_test(pi_gene, "nervous_system_process", "ZI", direction = "less")


###PI_win#################
#nerv sys dev
run_wilcoxon_test(pi_win, "nervous_system_development", "ZS", direction = "less")
run_wilcoxon_test(pi_win, "nervous_system_development", "ZH", direction = "less")
run_wilcoxon_test(pi_win, "nervous_system_development", "RAL", direction = "less")
run_wilcoxon_test(pi_win, "nervous_system_development", "FR", direction = "less")
run_wilcoxon_test(pi_win, "nervous_system_development", "ZI", direction = "less")

#sensory sys dev
run_wilcoxon_test(pi_win, "sensory_system_development", "ZS", direction = "less")
run_wilcoxon_test(pi_win, "sensory_system_development", "ZH", direction = "less")
run_wilcoxon_test(pi_win, "sensory_system_development", "RAL", direction = "less")
run_wilcoxon_test(pi_win, "sensory_system_development", "FR", direction = "less")
run_wilcoxon_test(pi_win, "sensory_system_development", "ZI", direction = "less")

#mechanosensory behavior
run_wilcoxon_test(pi_win, "mechanosensory_behavior", "ZS", direction = "less")
run_wilcoxon_test(pi_win, "mechanosensory_behavior", "ZH", direction = "less")
run_wilcoxon_test(pi_win, "mechanosensory_behavior", "RAL", direction = "less")
run_wilcoxon_test(pi_win, "mechanosensory_behavior", "FR", direction = "less")
run_wilcoxon_test(pi_win, "mechanosensory_behavior", "ZI", direction = "less")

#aggressive behavior
run_wilcoxon_test(pi_win, "aggressive_behavior", "ZS", direction = "less")
run_wilcoxon_test(pi_win, "aggressive_behavior", "ZH", direction = "less")
run_wilcoxon_test(pi_win, "aggressive_behavior", "RAL", direction = "less")
run_wilcoxon_test(pi_win, "aggressive_behavior", "FR", direction = "less")
run_wilcoxon_test(pi_win, "aggressive_behavior", "ZI", direction = "less")

#learning or memory
run_wilcoxon_test(pi_win, "learning_or_memory", "ZS", direction = "less")
run_wilcoxon_test(pi_win, "learning_or_memory", "ZH", direction = "less")
run_wilcoxon_test(pi_win, "learning_or_memory", "RAL", direction = "less")
run_wilcoxon_test(pi_win, "learning_or_memory", "FR", direction = "less")
run_wilcoxon_test(pi_win, "learning_or_memory", "ZI", direction = "less")

#nerv sys proc
run_wilcoxon_test(pi_win, "nervous_system_process", "ZS", direction = "less")
run_wilcoxon_test(pi_win, "nervous_system_process", "ZH", direction = "less")
run_wilcoxon_test(pi_win, "nervous_system_process", "RAL", direction = "less")
run_wilcoxon_test(pi_win, "nervous_system_process", "FR", direction = "less")
run_wilcoxon_test(pi_win, "nervous_system_process", "ZI", direction = "less")


###THETA#################
#nerv sys dev
run_wilcoxon_test(theta, "nervous_system_development", "ZS", direction = "less")
run_wilcoxon_test(theta, "nervous_system_development", "ZH", direction = "less")
run_wilcoxon_test(theta, "nervous_system_development", "RAL", direction = "less")
run_wilcoxon_test(theta, "nervous_system_development", "FR", direction = "less")
run_wilcoxon_test(theta, "nervous_system_development", "ZI", direction = "less")

#sensory sys dev
run_wilcoxon_test(theta, "sensory_system_development", "ZS", direction = "less")
run_wilcoxon_test(theta, "sensory_system_development", "ZH", direction = "less")
run_wilcoxon_test(theta, "sensory_system_development", "RAL", direction = "less")
run_wilcoxon_test(theta, "sensory_system_development", "FR", direction = "less")
run_wilcoxon_test(theta, "sensory_system_development", "ZI", direction = "less")

#mechanosensory behavior
run_wilcoxon_test(theta, "mechanosensory_behavior", "ZS", direction = "less")
run_wilcoxon_test(theta, "mechanosensory_behavior", "ZH", direction = "less")
run_wilcoxon_test(theta, "mechanosensory_behavior", "RAL", direction = "less")
run_wilcoxon_test(theta, "mechanosensory_behavior", "FR", direction = "less")
run_wilcoxon_test(theta, "mechanosensory_behavior", "ZI", direction = "less")

#aggressive behavior
run_wilcoxon_test(theta, "aggressive_behavior", "ZS", direction = "less")
run_wilcoxon_test(theta, "aggressive_behavior", "ZH", direction = "less")
run_wilcoxon_test(theta, "aggressive_behavior", "RAL", direction = "less")
run_wilcoxon_test(theta, "aggressive_behavior", "FR", direction = "less")
run_wilcoxon_test(theta, "aggressive_behavior", "ZI", direction = "less")

#learning or memory
run_wilcoxon_test(theta, "learning_or_memory", "ZS", direction = "less")
run_wilcoxon_test(theta, "learning_or_memory", "ZH", direction = "less")
run_wilcoxon_test(theta, "learning_or_memory", "RAL", direction = "less")
run_wilcoxon_test(theta, "learning_or_memory", "FR", direction = "less")
run_wilcoxon_test(theta, "learning_or_memory", "ZI", direction = "less")

#nerv sys proc
run_wilcoxon_test(theta, "nervous_system_process", "ZS", direction = "less")
run_wilcoxon_test(theta, "nervous_system_process", "ZH", direction = "less")
run_wilcoxon_test(theta, "nervous_system_process", "RAL", direction = "less")
run_wilcoxon_test(theta, "nervous_system_process", "FR", direction = "less")
run_wilcoxon_test(theta, "nervous_system_process", "ZI", direction = "less")


###TjD#################
#nerv sys dev
run_wilcoxon_test(tajimaD, "nervous_system_development", "ZS", direction = "less")
run_wilcoxon_test(tajimaD, "nervous_system_development", "ZH", direction = "less")
run_wilcoxon_test(tajimaD, "nervous_system_development", "RAL", direction = "less")
run_wilcoxon_test(tajimaD, "nervous_system_development", "FR", direction = "less")
run_wilcoxon_test(tajimaD, "nervous_system_development", "ZI", direction = "less")

#sensory sys dev
run_wilcoxon_test(tajimaD, "sensory_system_development", "ZS", direction = "less")
run_wilcoxon_test(tajimaD, "sensory_system_development", "ZH", direction = "less")
run_wilcoxon_test(tajimaD, "sensory_system_development", "RAL", direction = "less")
run_wilcoxon_test(tajimaD, "sensory_system_development", "FR", direction = "less")
run_wilcoxon_test(tajimaD, "sensory_system_development", "ZI", direction = "less")

#mechanosensory behavior
run_wilcoxon_test(tajimaD, "mechanosensory_behavior", "ZS", direction = "less")
run_wilcoxon_test(tajimaD, "mechanosensory_behavior", "ZH", direction = "less")
run_wilcoxon_test(tajimaD, "mechanosensory_behavior", "RAL", direction = "less")
run_wilcoxon_test(tajimaD, "mechanosensory_behavior", "FR", direction = "less")
run_wilcoxon_test(tajimaD, "mechanosensory_behavior", "ZI", direction = "less")

#aggressive behavior
run_wilcoxon_test(tajimaD, "aggressive_behavior", "ZS", direction = "less")
run_wilcoxon_test(tajimaD, "aggressive_behavior", "ZH", direction = "less")
run_wilcoxon_test(tajimaD, "aggressive_behavior", "RAL", direction = "less")
run_wilcoxon_test(tajimaD, "aggressive_behavior", "FR", direction = "less")
run_wilcoxon_test(tajimaD, "aggressive_behavior", "ZI", direction = "less")

#learning or memory
run_wilcoxon_test(tajimaD, "learning_or_memory", "ZS", direction = "less")
run_wilcoxon_test(tajimaD, "learning_or_memory", "ZH", direction = "less")
run_wilcoxon_test(tajimaD, "learning_or_memory", "RAL", direction = "less")
run_wilcoxon_test(tajimaD, "learning_or_memory", "FR", direction = "less")
run_wilcoxon_test(tajimaD, "learning_or_memory", "ZI", direction = "less")

#nerv sys proc
run_wilcoxon_test(tajimaD, "nervous_system_process", "ZS", direction = "less")
run_wilcoxon_test(tajimaD, "nervous_system_process", "ZH", direction = "less")
run_wilcoxon_test(tajimaD, "nervous_system_process", "RAL", direction = "less")
run_wilcoxon_test(tajimaD, "nervous_system_process", "FR", direction = "less")
run_wilcoxon_test(tajimaD, "nervous_system_process", "ZI", direction = "less")



###Plotting########################################################################################################################

plot_full_data <- function(data, go_terms, value_columns, dataset_name, statistic, test_group_color, sig_labels = NULL) {
  
  all_rows <- list()  # empty list to store all the reshaped chunks
  
  for (pop in value_columns) {
    for (go in go_terms) {
      
      go_label <- str_replace_all(go, "_", " ")
      pop_values <- data[[pop]]
      go_membership <- data[[go]]
      
      df_chunk <- tibble(
        Value = pop_values,
        Population = pop,
        GO = go_label,
        Group = if_else(go_membership == go_label, "In GO", "Not in GO")
      )
      
      all_rows[[length(all_rows) + 1]] <- df_chunk
    }
  }
  
  # Combine all chunks into one long tidy dataframe
  plot_data <- bind_rows(all_rows)
  
  # Set the facet order
  plot_data$Population <- factor(plot_data$Population, levels = value_columns)
  plot_data$GO <- factor(plot_data$GO, levels = str_replace_all(go_terms, "_", " "))

  
  # Plot
  p <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group)) +
    #stat_summary(fun = mean, geom = "bar", position = position_dodge(0.9), width = 0.4) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5, position = position_dodge(0.5)) +
    stat_summary(fun.data = function(x) mean_sdl(x, mult = 1), geom = "errorbar", position = position_dodge(0.9), width = 0.2, color = "black") +
    facet_grid(rows = vars(GO), cols = vars(Population), scales = "free_y") +
    scale_fill_manual(values = c("In GO" = test_group_color, "Not in GO" = "gray70")) +
    theme_bw(base_size = 10) +
    labs(title = dataset_name, x = NULL, y = statistic, fill = "Group") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(0.5, "lines"),
      plot.margin = margin(5, 5, 5, 10),
      legend.position = "right"
    )
  
  # Add significance annotations
  if (!is.null(sig_labels)) {
    sig_labels_annot <- sig_labels %>%
      mutate(GO = str_replace_all(GO, "_", " "))  # match the facet labels
    
    # Estimate y-position for annotation
    summary_df <- plot_data %>%
      group_by(GO, Population, Group) %>%
      summarise(
        mean = mean(Value, na.rm = TRUE),
        sd = sd(Value, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(top = mean + sd)
    
    # Get row-wise (GO-wise) max
    row_max <- summary_df %>%
      group_by(GO) %>%
      summarise(y_pos = max(top), .groups = "drop")
    
    # Join to sig_labels
    label_df <- sig_labels_annot %>%
      mutate(GO = str_replace_all(GO, "_", " ")) %>%
      inner_join(row_max, by = "GO")
    
    label_df$GO <- factor(label_df$GO, levels = levels(plot_data$GO))
    label_df$Population <- factor(label_df$Population, levels = levels(plot_data$Population))
    
    p <- p + geom_text(data = label_df, 
                      aes(x = 1.5, y = y_pos, label = label), 
                      inherit.aes = FALSE, size = 8, fontface = "bold") +
                      scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
                      coord_cartesian(clip = "off") +
                      theme(plot.margin = margin(15, 10, 10, 10))
  }
  
  p

}


go_terms <- c(
  "nervous_system_development",
  "sensory_system_development",
  "mechanosensory_behavior",
  "aggressive_behavior",
  "learning_or_memory",
  "nervous_system_process"
)
pop_comp_cols <- c("ZS.vs.RAL", "ZS.vs.FR", "ZS.vs.ZI", "ZH.vs.RAL", "ZH.vs.FR", "ZH.vs.ZI", "RAL.vs.FR", "FR.vs.ZI", "ZI.vs.RAL")
pop_cols <- c("ZS", "ZH", "RAL", "FR", "ZI")



fst_gene_sig <- tibble(
  GO = rep("nervous_system_process", 7),
  Population = c("ZS.vs.RAL", "ZS.vs.FR", "ZH.vs.RAL", "ZH.vs.FR", "ZH.vs.ZI", "FR.vs.ZI", "ZI.vs.RAL"),
  label = c("***", "**", "**", "**", "*", "**", "**")
)


fst_win_sig <- tibble(
  GO = c(rep("nervous_system_development", 2), rep("sensory_system_development", 2), "aggressive_behavior", rep("learning_or_memory", 8), rep("nervous_system_process", 8)),
  Population = c("ZS.vs.ZI", "ZH.vs.RAL", "ZH.vs.RAL", "ZH.vs.FR", "ZS.vs.RAL", rep(c("ZS.vs.RAL", "ZS.vs.FR", "ZS.vs.ZI", "ZH.vs.RAL", "ZH.vs.FR", "ZH.vs.ZI", "FR.vs.ZI", "ZI.vs.RAL"), 2)),
  label = c("*", "*", "*", "*", "*", rep("***", 4), rep(c("**", "***"), 2), rep("***", 4), "*", rep("***", 3))
)


pi_gene_sig <- tibble(
  GO = c(rep("nervous_system_development", 5), rep("sensory_system_development", 5), rep("learning_or_memory", 5)),
  Population = rep(c("ZS", "ZH", "RAL", "FR", "ZI"), 3),
  label = rep("***", 15)
)


pi_win_sig <- tibble(
  GO = c(rep("mechanosensory_behavior", 5), rep("learning_or_memory", 2), "nervous_system_process"),
  Population = c("ZS", "ZH", "RAL", "FR", "ZI", "RAL", "FR", "RAL"),
  label = c(rep("***", 7), "**")
)


theta_sig <- tibble(
  GO = c(rep("mechanosensory_behavior", 5), rep("learning_or_memory", 2), "nervous_system_process"),
  Population = c("ZS", "ZH", "RAL", "FR", "ZI", "RAL", "FR", "RAL"),
  label = c(rep("***", 7), "**")
)


tjd_sig <- tibble(
  GO = c("nervous_system_development", rep("sensory_system_development", 2), rep("mechanosensory_behavior", 4), rep("learning_or_memory", 3), "nervous_system_process"),
  Population = c("ZI", "ZH", "ZI", "ZS", "ZH", "RAL", "FR", "ZH", "RAL", "ZI", "ZI"),
  label = c(rep("**", 2), rep("***", 6), "*", "***", "**")
)


plot_full_data(fst_gene, go_terms, pop_comp_cols, "FST (Gene-level)", "FST", "lightseagreen")
plot_full_data(fst_win, go_terms, pop_comp_cols, "FST (10kbp Windows)", "FST", "lightseagreen")
plot_full_data(pi_gene, go_terms, pop_cols, "Pi (Gene-level)", "Pi", "coral", sig_labels = pi_gene_sig)
plot_full_data(pi_win, go_terms, pop_cols, "Pi (10kbp Windows)", "Pi", "coral", sig_labels = pi_win_sig)
plot_full_data(theta, go_terms, pop_cols, "Watterson's Theta (10kbp Windows)", "Theta", "firebrick1")
plot_full_data(tajimaD, go_terms, pop_cols, "Tajima's D (10kbp Windows)", "TjD", "darkmagenta")















plot_full_data_highlighted <- function(data, go_terms, value_columns, dataset_name, statistic, test_group_color, background_color, sig_labels) {

  
  # Build long-format plot_data
  all_rows <- list()
  for (pop in value_columns) {
    for (go in go_terms) {
      go_label <- str_replace_all(go, "_", " ")
      df_chunk <- tibble(
        Value = data[[pop]],
        Population = pop,
        GO = go_label,
        Group = if_else(data[[go]] == go_label, "In GO", "Not in GO")
      )
      all_rows[[length(all_rows) + 1]] <- df_chunk
    }
  }
  
  plot_data <- bind_rows(all_rows)
  
  # Set order
  plot_data$Population <- factor(plot_data$Population, levels = value_columns)
  plot_data$GO <- factor(plot_data$GO, levels = str_replace_all(go_terms, "_", " "))
  
  # Format significance label dataframe
  sig_labels_annot <- sig_labels %>%
    mutate(GO = str_replace_all(GO, "_", " "))
  
  sig_labels_annot$Population <- factor(sig_labels_annot$Population, levels = levels(plot_data$Population))
  sig_labels_annot$GO <- factor(sig_labels_annot$GO, levels = levels(plot_data$GO))
  
  # Build highlight rectangles
  highlight_df <- expand.grid(
    xmin = 0.5,
    xmax = 2.5,
    ymin = -Inf,
    ymax = Inf
  ) %>%
    crossing(sig_labels_annot)
  
  # Calculate per-row max for star placement
  summary_df <- plot_data %>%
    group_by(GO, Population, Group) %>%
    summarise(mean = mean(Value, na.rm = TRUE),
              sd = sd(Value, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(top = mean + sd)
  
  row_max <- summary_df %>%
    group_by(GO) %>%
    summarise(y_pos = max(top), .groups = "drop")
  
  label_df <- sig_labels_annot %>%
    inner_join(row_max, by = "GO")
  
  # Plot
  go_display_labels <- c(
    "nervous system development" = "nervous sys. dev",
    "sensory system development" = "sensory sys. dev",
    "mechanosensory behavior" = "mechanosensory behavior",
    "aggressive behavior" = "aggressive behavior",
    "learning or memory" = "learning or memory",
    "nervous system process" = "nervous sys. process"
  )
  
  p <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group)) +
    
    # Background first
    geom_rect(
      data = highlight_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = background_color,
      alpha = 0.2
    ) +
    
    # Main data
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3.5, position = position_dodge(0.5)) +
    stat_summary(fun.data = function(x) mean_sdl(x, mult = 1), geom = "errorbar", 
                 width = 0.2, color = "black", position = position_dodge(0.9)) +
    
    # Stars
    geom_text(data = label_df,
              aes(x = 1.5, y = y_pos, label = label),
              inherit.aes = FALSE,
              size = 6, fontface = "bold") +
    
    facet_grid(rows = vars(GO), cols = vars(Population), scales = "free_y", labeller = labeller(GO = go_display_labels)) +
    scale_fill_manual(values = c("In GO" = test_group_color, "Not in GO" = "gray70")) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(title = dataset_name, x = NULL, y = statistic, fill = "Group") +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(0.5, "lines"),
      plot.margin = margin(15, 10, 10, 10),
      legend.position = "right"
    )
  
  return(p)
}




# FST (Gene-level)
plot_full_data_highlighted(
  fst_gene, go_terms, pop_comp_cols,
  dataset_name = "FST (Gene-level)",
  statistic = "FST",
  test_group_color = "lightseagreen",
  background_color = "greenyellow",
  sig_labels = fst_gene_sig
)

# FST (10kbp Windows)
plot_full_data_highlighted(
  fst_win, go_terms, pop_comp_cols,
  dataset_name = "FST (10kbp Windows)",
  statistic = "FST",
  test_group_color = "lightseagreen",
  background_color = "greenyellow",
  sig_labels = fst_win_sig
)

# Pi (Gene-level)
plot_full_data_highlighted(
  pi_gene, go_terms, pop_cols,
  dataset_name = "Pi (Gene-level)",
  statistic = "Pi",
  test_group_color = "coral",
  background_color = "lightpink",
  sig_labels = pi_gene_sig
)

# Pi (10kbp Windows)
plot_full_data_highlighted(
  pi_win, go_terms, pop_cols,
  dataset_name = "Pi (10kbp Windows)",
  statistic = "Pi",
  test_group_color = "coral",
  background_color = "lightpink",
  sig_labels = pi_win_sig
)

# Watterson's Theta
plot_full_data_highlighted(
  theta, go_terms, pop_cols,
  dataset_name = "Watterson's Theta (10kbp Windows)",
  statistic = "Theta",
  test_group_color = "firebrick1",
  background_color = "lightpink",
  sig_labels = theta_sig
)

# Tajima's D
plot_full_data_highlighted(
  tajimaD, go_terms, pop_cols,
  dataset_name = "Tajima's D (10kbp Windows)",
  statistic = "Tajima's D",
  test_group_color = "darkmagenta",
  background_color = "lightpink",
  sig_labels = tjd_sig
)


