#!/bin/bash

vcftools --vcf ZS_biallelicSNPs_dmel_r6_allChr.vcf --SNPdensity 10000 --out ZS
vcftools --vcf ZH_biallelicSNPs_dmel_r6_allChr.vcf --SNPdensity 10000 --out ZH
vcftools --vcf RAL_biallelicSNPs_dmel_r6_allChr.vcf --SNPdensity 10000 --out RAL
vcftools --vcf FR_biallelicSNPs_dmel_r6_allChr.vcf --SNPdensity 10000 --out FR
vcftools --vcf ZI_biallelicSNPs_dmel_r6_allChr.vcf --SNPdensity 10000 --out ZI
