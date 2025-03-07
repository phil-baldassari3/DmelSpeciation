import pandas as pd

"""
#open SRA csv
df = pd.read_csv("sample_SRAs.csv")

#list of SRAs
SRAs = df["SRA"].to_list()
"""

with open("SRAs_to_download.txt", "r") as sf:
    SRAstr = sf.read()
    SRAs = SRAstr.split("\n")


#empty commands string
commands_str = ""

#looping through SRAs
for sra in SRAs:
    accession = sra

    #commands_str += "(prefetch {acc}; fasterq-dump --split-3 --skip-technical {acc}) && (rm -r {acc}) && (gzip {acc}_1.fastq & gzip {acc}_2.fastq) && (mv {acc}_1.fastq.gz /Volumes/External\ Files/dmel_SRAs) && (mv {acc}_2.fastq.gz /Volumes/External\ Files/dmel_SRAs)\n".format(acc=accession)
    commands_str += "(prefetch {acc}; fasterq-dump --split-3 --skip-technical {acc}) && (rm -r {acc}) && (gzip {acc}_1.fastq & gzip {acc}_2.fastq)\n".format(acc=accession)


with open("download_SRAs_to_owlsnest.sh", "w") as commands:
    commands.write(commands_str)