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

    commands_str += "java -jar $PICARD MarkDuplicates INPUT=sorted_{name}.sam OUTPUT=deduped/deduped_{name}.sam METRICS_FILE=deduped/metrics/{name}_metrics.txt\n".format(name=SRA2name[accession])


with open("MarkDuplicates_commands.txt", "w") as commands:
    commands.write(commands_str)