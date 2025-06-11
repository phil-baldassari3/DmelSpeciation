#!/bin/bash

vcftools --vcf ZS_dmel_r6_allChr.vcf --positions include_sites_ZS.txt --recode --out ZS_biallelicSNPs_dmel_r6_allChr
vcftools --vcf ZH_dmel_r6_allChr.vcf --positions include_sites_ZH.txt --recode --out ZH_biallelicSNPs_dmel_r6_allChr
vcftools --vcf RAL_dmel_r6_allChr.vcf --positions include_sites_RAL.txt --recode --out RAL_biallelicSNPs_dmel_r6_allChr
vcftools --vcf FR_dmel_r6_allChr.vcf --positions include_sites_FR.txt --recode --out FR_biallelicSNPs_dmel_r6_allChr
vcftools --vcf ZI_dmel_r6_allChr.vcf --positions include_sites_ZI.txt --recode --out ZI_biallelicSNPs_dmel_r6_allChr
