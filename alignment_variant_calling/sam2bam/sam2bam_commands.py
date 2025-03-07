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

    commands_str += 'java -jar $PICARD SamFormatConverter INPUT={name}.sam OUTPUT=sam2bam/{name}.bam\n'.format(name=SRA2name[accession])


with open("sam2bam_commands.txt", "w") as commands:
    commands.write(commands_str)