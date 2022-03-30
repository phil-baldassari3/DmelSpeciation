#!/bin/bash

# make fasta alignment files
cd /Users/philipbaldassari/Desktop/zim-cos_ChrXseq

python seq2fas_ChrX.py

cd ../zim-cos_Chr2Lseq

python seq2fas_Chr2L.py

cd ../zim-cos_Chr2Rseq

python seq2fas_Chr2R.py

cd ../zim-cos_Chr3Lseq

python seq2fas_Chr3L.py

cd ../zim-cos_Chr3Rseq

python seq2fas_Chr3R.py

echo "Individual fasta files created"

# Make vcfs

# create ChrX vcf
cd ../zim-cos_ChrX

cat $(ls) > zim-cos_ChrX.fa

for fas in *.fas
do
	rm $fas
done

snp-sites -v -o pre_zim-cos_ChrX.vcf zim-cos_ChrX.fa

rm zim-cos_ChrX.fa

echo "1 X" >> chrX_name_conv.txt
bcftools annotate --rename-chrs chrX_name_conv.txt pre_zim-cos_ChrX.vcf > ../zim-cos_vcfs/chr_pre_zim-cos_ChrX.vcf

rm pre_zim-cos_ChrX.vcf

echo "vcf for Chrom X made"

# create Chr2L vcf
cd ../zim-cos_Chr2L

cat $(ls) > zim-cos_Chr2L.fa

for fas in *.fas
do
	rm $fas
done

snp-sites -v -o pre_zim-cos_Chr2L.vcf zim-cos_Chr2L.fa

rm zim-cos_Chr2L.fa

echo "1 2L" >> chr2L_name_conv.txt
bcftools annotate --rename-chrs chr2L_name_conv.txt pre_zim-cos_Chr2L.vcf > ../zim-cos_vcfs/chr_pre_zim-cos_Chr2L.vcf

rm pre_zim-cos_Chr2L.vcf

echo "vcf for Chrom 2L made"

# create Chr2R vcf
cd ../zim-cos_Chr2R

cat $(ls) > zim-cos_Chr2R.fa

for fas in *.fas
do
	rm $fas
done

snp-sites -v -o pre_zim-cos_Chr2R.vcf zim-cos_Chr2R.fa

rm zim-cos_Chr2R.fa

echo "1 2R" >> chr2R_name_conv.txt
bcftools annotate --rename-chrs chr2R_name_conv.txt pre_zim-cos_Chr2R.vcf > ../zim-cos_vcfs/chr_pre_zim-cos_Chr2R.vcf

rm pre_zim-cos_Chr2R.vcf

echo "vcf for Chrom 2R made"

# create Chr3L vcf
cd ../zim-cos_Chr3L

cat $(ls) > zim-cos_Chr3L.fa

for fas in *.fas
do
	rm $fas
done

snp-sites -v -o pre_zim-cos_Chr3L.vcf zim-cos_Chr3L.fa

rm zim-cos_Chr3L.fa

echo "1 3L" >> chr3L_name_conv.txt
bcftools annotate --rename-chrs chr3L_name_conv.txt pre_zim-cos_Chr3L.vcf > ../zim-cos_vcfs/chr_pre_zim-cos_Chr3L.vcf

rm pre_zim-cos_Chr3L.vcf

echo "vcf for Chrom 3L made"

# create Chr3R vcf
cd ../zim-cos_Chr3R

cat $(ls) > zim-cos_Chr3R.fa

for fas in *.fas
do
	rm $fas
done

snp-sites -v -o pre_zim-cos_Chr3R.vcf zim-cos_Chr3R.fa

rm zim-cos_Chr3R.fa

echo "1 3R" >> chr3R_name_conv.txt
bcftools annotate --rename-chrs chr3R_name_conv.txt pre_zim-cos_Chr3R.vcf > ../zim-cos_vcfs/chr_pre_zim-cos_Chr3R.vcf

rm pre_zim-cos_Chr3R.vcf

echo "vcf for Chrom 3R made"


# changing directory to vcf directory
cd ../zim-cos_vcfs

#make diploid
python hap2dip.py
rm chr_pre_zim-cos_ChrX.vcf chr_pre_zim-cos_Chr2L.vcf chr_pre_zim-cos_Chr2R.vcf chr_pre_zim-cos_Chr3L.vcf chr_pre_zim-cos_Chr3R.vcf

echo "Made vcfs diploid (assuming inbred lines are almost entirely nomozygous and heterozygous regions have been masked)"
echo "temporary files removed"

# remove sites w/ N as common ALT
for file in *.vcf
do
	grep -v "*," $file > N_uncommon/N_uncommon_$file
done

echo "made files with N as uncommon ALT"

#keep sites w/ N as common ALT, header is also removed
for f in *.vcf
do
	grep "*," $f > N_common/N_common_$f
done

echo "made files with N as common ALT"

# Remove temporary files
rm *.vcf

echo "removed temprorary files"

# run vcf fixer scripts
cd N_common
python vcf_fixer_N_common.py
mv fixed_*.vcf ..
rm *.vcf

cd ../N_uncommon
python vcf_fixer_N_uncommon.py
mv fixed_*.vcf ..
rm *.vcf

echo "made fixed vcf that need to be concatenated"

cd ..

# concatenating vcfs
cat fixed_N_uncommon_dip_chr_pre_zim-cos_ChrX.vcf fixed_N_common_dip_chr_pre_zim-cos_ChrX.vcf > unsorted_zim-cos_ChrX.vcf
cat fixed_N_uncommon_dip_chr_pre_zim-cos_Chr2L.vcf fixed_N_common_dip_chr_pre_zim-cos_Chr2L.vcf > unsorted_zim-cos_Chr2L.vcf
cat fixed_N_uncommon_dip_chr_pre_zim-cos_Chr2R.vcf fixed_N_common_dip_chr_pre_zim-cos_Chr2R.vcf > unsorted_zim-cos_Chr2R.vcf
cat fixed_N_uncommon_dip_chr_pre_zim-cos_Chr3L.vcf fixed_N_common_dip_chr_pre_zim-cos_Chr3L.vcf > unsorted_zim-cos_Chr3L.vcf
cat fixed_N_uncommon_dip_chr_pre_zim-cos_Chr3R.vcf fixed_N_common_dip_chr_pre_zim-cos_Chr3R.vcf > unsorted_zim-cos_Chr3R.vcf

rm fixed_*.vcf

echo "vcf concatnetated, temporary files removed"

# sort vcfs
bcftools sort unsorted_zim-cos_ChrX.vcf > zim-cos_ChrX.vcf
bcftools sort unsorted_zim-cos_Chr2L.vcf > zim-cos_Chr2L.vcf
bcftools sort unsorted_zim-cos_Chr2R.vcf > zim-cos_Chr2R.vcf
bcftools sort unsorted_zim-cos_Chr3L.vcf > zim-cos_Chr3L.vcf
bcftools sort unsorted_zim-cos_Chr3R.vcf > zim-cos_Chr3R.vcf

rm unsorted_*.vcf

echo "vcfs sorted, temporary files removed"

echo "DONE!"












