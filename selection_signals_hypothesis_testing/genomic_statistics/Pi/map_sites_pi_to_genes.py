import pandas as pd

#open data tables
pi = pd.read_csv("selected_pops_sites_pi.csv")
gene_map = pd.read_csv("biallelicSNPs_missing0.05_dmel_r6_allChr_gene_ann.tsv", sep="\t", keep_default_na=False)

#these effects demostrate the variant is not inside the gene
excluded_effects = ["intergenic_region", "upstream_gene_variant", "downstream_gene_variant"]

#filter gene map
filtered_gene_map = gene_map[(gene_map["ANN[*].BIOTYPE"] == "protein_coding") & ~gene_map["ANN[*].EFFECT"].isin(excluded_effects)]
filtered_gene_map = filtered_gene_map[["CHROM", "POS", "ANN[*].GENEID", "ANN[*].GENE"]]
filtered_gene_map = filtered_gene_map.drop_duplicates()


#making gene table
gene_table = filtered_gene_map.rename(columns={"ANN[*].GENEID": "FBgn", "ANN[*].GENE": "gene_symbol"})
gene_table = gene_table.groupby(["CHROM", "POS"], as_index=False).agg({"FBgn": lambda x: "; ".join(x), "gene_symbol": lambda x: "; ".join(x)})


#merge FST and gene table
merged = pi.merge(gene_table, on=["CHROM", "POS"], how="left")
merged["FBgn"] = merged["FBgn"].fillna("None")
merged["gene_symbol"] = merged["gene_symbol"].fillna("None")

#saving mapped table
merged.to_csv("selected_pops_sites_pi_mapped2genes.csv", index=False)