#!/bin/bash

vcftools --vcf ZS_NoSingletons_dmel_r6_allChr.vcf --freq --out ZS
vcftools --vcf ZH_NoSingletons_dmel_r6_allChr.vcf --freq --out ZH
vcftools --vcf RAL_NoSingletons_dmel_r6_allChr.vcf --freq --out RAL
vcftools --vcf FR_NoSingletons_dmel_r6_allChr.vcf --freq --out FR
vcftools --vcf ZI_NoSingletons_dmel_r6_allChr.vcf --freq --out ZI