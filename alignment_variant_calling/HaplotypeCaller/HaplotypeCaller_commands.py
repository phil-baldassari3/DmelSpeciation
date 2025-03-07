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

    commands_str += 'gatk --java-options "-Xmx8g" HaplotypeCaller -R /home/tuk40537/scratch/Dmel_reads/reference_genome_gatk/dmel-all-chromosome-r6.56.fasta -I {name}.bam -O HaplotypeCaller/{name}.g.vcf.gz -ERC GVCF --tmp-dir /local_scratch/tmp/$PBS_JOBID\n'.format(name=SRA2name[accession])


with open("HaplotypeCaller_commands.txt", "w") as commands:
    commands.write(commands_str)