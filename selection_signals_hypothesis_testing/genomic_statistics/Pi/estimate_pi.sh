#!/bin/bash

vcftools --vcf ZS_dmel_r6_allChr.vcf --site-pi --out ZS
vcftools --vcf ZH_dmel_r6_allChr.vcf --site-pi --out ZH
vcftools --vcf RAL_dmel_r6_allChr.vcf --site-pi --out RAL
vcftools --vcf FR_dmel_r6_allChr.vcf --site-pi --out FR
vcftools --vcf ZI_dmel_r6_allChr.vcf --site-pi --out ZI


vcftools --vcf ZS_dmel_r6_allChr.vcf --window-pi 10000 --out ZS
vcftools --vcf ZH_dmel_r6_allChr.vcf --window-pi 10000 --out ZH
vcftools --vcf RAL_dmel_r6_allChr.vcf --window-pi 10000 --out RAL
vcftools --vcf FR_dmel_r6_allChr.vcf --window-pi 10000 --out FR
vcftools --vcf ZI_dmel_r6_allChr.vcf --window-pi 10000 --out ZI
