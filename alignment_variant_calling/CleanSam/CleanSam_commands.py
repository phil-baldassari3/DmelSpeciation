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

    commands_str += "java -jar $PICARD CleanSam INPUT=RGadded_{name}.sam OUTPUT=cleaned/cleaned_{name}.sam\n".format(name=SRA2name[accession])


with open("CleanSam_commands.txt", "w") as commands:
    commands.write(commands_str)