#!/bin/bash

vcftools --vcf ZS_dmel_r6_allChr.vcf --freq --out ZS_freq4filtering
vcftools --vcf ZH_dmel_r6_allChr.vcf --freq --out ZH_freq4filtering
vcftools --vcf RAL_dmel_r6_allChr.vcf --freq --out RAL_freq4filtering
vcftools --vcf FR_dmel_r6_allChr.vcf --freq --out FR_freq4filtering
vcftools --vcf ZI_dmel_r6_allChr.vcf --freq --out ZI_freq4filtering
