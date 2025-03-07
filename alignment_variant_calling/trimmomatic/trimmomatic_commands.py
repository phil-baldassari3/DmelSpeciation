import pandas as pd


#open SRA csv
df = pd.read_csv("sample_SRAs.csv")

#list of SRAs
SRAscsv = df["SRA"].to_list()
#list of names
Names = df["Name"].to_list()

with open("SRAs.txt", "r") as sf:
    SRAstr = sf.read()
    SRAs = SRAstr.split("\n")


#SRA to Name dict
SRA2name = dict(zip(SRAscsv, Names))

#empty commands string
commands_str = ""

#looping through SRAs
for sra in SRAs:
    accession = sra

    commands_str += "java -jar /home/tuk40537/Trimmomatic-0.39/trimmomatic-0.39.jar PE {acc}_1.fastq {acc}_2.fastq trimmed/{name}_1P.fastq singletons/{name}_1U.fastq trimmed/{name}_2P.fastq singletons/{name}_2U.fastq ILLUMINACLIP:/home/tuk40537/Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\n".format(acc=accession, name=SRA2name[accession])


with open("trimmomatic_commands.txt", "w") as commands:
    commands.write(commands_str)