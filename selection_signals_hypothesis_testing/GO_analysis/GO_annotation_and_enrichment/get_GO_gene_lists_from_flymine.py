from intermine.webservice import Service
import pandas as pd


service = Service("https://www.flymine.org/flymine/service")

def flymineGO2df(Goterm):

    query = service.new_query("Gene")

    query.add_constraint("goAnnotation.ontologyTerm.parents", "GOTerm")
    query.add_constraint("goAnnotation.ontologyTerm", "GOTerm")

    query.add_view(
        "primaryIdentifier", "symbol", "goAnnotation.ontologyTerm.parents.name",
        "goAnnotation.ontologyTerm.parents.identifier"
    )

    query.add_sort_order("Gene.secondaryIdentifier", "ASC")

    query.add_constraint("goAnnotation.ontologyTerm.parents.name", "=", Goterm, code="A")
    query.add_constraint("organism.name", "=", "Drosophila melanogaster", code="B")

    #getting df
    df = query.dataframe()

    #extracting columns
    df = df[["Gene.primaryIdentifier", "Gene.symbol", "Gene.goAnnotation.ontologyTerm.parents.name", "Gene.goAnnotation.ontologyTerm.parents.identifier"]]

    return df



def gene_list_from_GOdf(df):

    genes = df["Gene.primaryIdentifier"].to_list()
    genes = list(set(genes))

    return genes




"""
listofGOterms = [
    "sensory system development", "digestive system development",
    "renal system development", "urogenital system development", "circulatory system development",
    "exocrine system development", "respiratory system development",
    "nervous system development", "reproductive system development", "immune system development"]

listofGOterms = ["locomotory behavior", 
"negative regulation of behavior", "thermosensory behavior",
"adult behavior", "aggressive behavior", 
"rhythmic behavior", 
"larval behavior", "regulation of behavior", "mechanosensory behavior", 
"positive regulation of behavior",
"reproductive behavior",
"visual behavior", "feeding behavior", 
"chemosensory behavior", "learning or memory"
]


listofGOterms = ["copulation", "female mating behavior",
"male mating behavior", "courtship behavior"
]
"""

listofGOterms = [
    "excretion", "hepaticobiliary system process", "muscle system process",
    "renal system process", "circulatory system process", "endocrine process",
    "respiratory system process", "digestive system process", "nervous system process"
]


for go in listofGOterms:

    GOdf = flymineGO2df(go)
    GOgenes = gene_list_from_GOdf(GOdf)

    name = go.replace(" ", "_")

    with open(f"{name}_genes.txt", "w") as genefile:
        genefile.write("\n".join(GOgenes))
