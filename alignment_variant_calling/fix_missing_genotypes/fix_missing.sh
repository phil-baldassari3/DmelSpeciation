#!/bin/sh

bcftools +setGT $1 -Oz -o FixMissing/temp_$1 - -- -t q -n . -i 'FORMAT/DP=0 | SMPL_MAX(FORMAT/PL)=0'

zcat FixMissing/temp_$1 | sed -E 's/\.\/\.\:[^\t]*:[^\t]*:[^\t]*:[^\t]*/\.\/\.\:\.\:\.\:\.\:\./g' | gzip > FixMissing/fixed_$1