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

    commands_str += "java -jar $PICARD AddOrReplaceReadGroups INPUT={name}.sam OUTPUT=AddReadGroups/RGadded_{name}.sam RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM={name}\n".format(name=SRA2name[accession])


with open("AddReadGroups_commands.txt", "w") as commands:
    commands.write(commands_str)