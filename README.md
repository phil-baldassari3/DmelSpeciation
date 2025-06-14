# DmelSpeciation

Scripts for the *D. melanogaster* Zimbabwe speciation project.

Directories and Scripts:

- **loci_table_toolkit**: Contains a python module and its yml env file for annotating Dmel genomic data with genes and GO terms and analyzing enrichment

- **alignment_variant_calling**: BWA -> PicardTools -> GATK Pipeline for aligning and variant calling Dmel SRA reads

- **pop_structure_analysis**: Scripts for PCA and ADMIXTURE analysis
    - `samples_map.R`: plot sample populations on a map
    - `pca-plot.R`: PCA of Dmel samples
    - `PCA_kmeans.py`: cluster PCA with kmeans clustering
    - `pca-kmeans.R`: plot PCA with kmeans clusters
    - `admixture_commands.txt`: run ADMIXTURE
    - `plot_ADMIXTURE.R`: plot structure plot

- **FST_estimation**: Estimating FST and various analyses
    - **estimating_Hudson_FST**: Scripts to estimate Hudson FST
    - **sliding_window_FST**: Computing FST in sliding windows
    - **Zim_specific_SNPs**: Categorizing Zim-specific SNPs
    - **speciation_islands**: Finding "Islands of Speciation"
    - **fst4pg_analysis**: Running fst4pg to find highly differentiated genomic regions
    - `full_circos.R`: plotting Circos plots

- **Hudson_Fst_from_vcftools_frq**: Estimating Hudson FST from allele frequency vcftools output
    - `merge_frq_files.py`: merge vcftools allele frequency outputs
    - `estimate_Hudson_FST.py`: compute Hudson FST from frq files

- **find_candidate_genes**: Scripts for estimating per gene FST and finding candidate genes
    - `map2genes.py`: map FST data to genes
    - `avgFst_by_gene_and_comparison.py`: compute average FST per gene

- **selection_signals_hypothesis_testing**: Testing hypothesis of finding selection signals in specific populatiosn and GO categories
    - **preprocessing_vcfs**: Scripts for filtering raw vcfools and subsetting selected populations
    - **genomic_statistics**: estimating FST, Pi, Watterson's Theta, and Tajima's D and mapping to genes
    - **GO_analysis**: Annotating genes with GO terms and running GO enrichment analysis
    - `runPCAdapt.R`: 
