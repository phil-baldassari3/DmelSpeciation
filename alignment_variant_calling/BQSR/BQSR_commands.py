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

    commands_str += 'gatk --java-options "-Xmx8g" BaseRecalibrator -I fixedmates_{name}.sam -R /home/tuk40537/scratch/Dmel_reads/reference_genome_gatk/dmel-all-chromosome-r6.56.fasta --known-sites /home/tuk40537/scratch/Dmel_reads/known_variants/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.vcf.gz --tmp-dir /local_scratch/tmp/$PBS_JOBID -O tables4BQSR/{name}.table\n'.format(name=SRA2name[accession])


with open("BQSR_commands.txt", "w") as commands:
    commands.write(commands_str)