from intermine.webservice import Service
import os

service = Service("https://www.flymine.org/flymine/service", token="H109a6MbtcXbrct6wczb")


def upload_gene_list(directory):
    """
    Function takes in a directory path of gene list files (.txt one gene per line)
    and uploads the list to Flymine using the Intermine API. Set service and token above.
    """

    for file in os.listdir(directory):
        if file.endswith(".txt"):
            genefile = os.path.join(directory, file)
            listname = file.replace(".txt", "")

            service.create_list(genefile, "Gene", listname, "Drosophila melanogaster")




upload_gene_list("fst_top10/SNPavg_fst")
upload_gene_list("fst_top10/pergene_fst")
upload_gene_list("fst_top10/windowed_fst")
upload_gene_list("pi_bottom10/pergene_pi")
upload_gene_list("pi_bottom10/windowed_pi")
upload_gene_list("theta_bottom10")
upload_gene_list("tajimasD_bottom10")

upload_gene_list("fst_top5/SNPavg_fst")
upload_gene_list("fst_top5/pergene_fst")
upload_gene_list("fst_top5/windowed_fst")
upload_gene_list("pi_bottom5/pergene_pi")
upload_gene_list("pi_bottom5/windowed_pi")
upload_gene_list("theta_bottom5")
upload_gene_list("tajimasD_bottom5")

upload_gene_list("fst_top1/SNPavg_fst")
upload_gene_list("fst_top1/pergene_fst")
upload_gene_list("fst_top1/windowed_fst")
upload_gene_list("pi_bottom1/pergene_pi")
upload_gene_list("pi_bottom1/windowed_pi")
upload_gene_list("theta_bottom1")
upload_gene_list("tajimasD_bottom1")