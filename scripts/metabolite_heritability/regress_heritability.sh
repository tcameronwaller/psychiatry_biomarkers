#!/bin/bash

#chmod u+x script.sh

# Read private, local file paths.
echo "read private file path variables..."
cd ~/paths
path_temporary=$(<"./temporary_bipolar_metabolism.txt")
path_waller="$path_temporary/waller"
path_access="$path_waller/dock/access"
path_disequilibrium="$path_access/disequilibrium"
path_alleles="$path_access/alleles"
path_metabolites="$path_access/metabolites"
path_heritability="$path_waller/dock/heritability"
path_ldsc=$(<"./tools_ldsc.txt")

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

# M00053 glutamine
# -- 1030131 SNPs
# -- heritability 0.4114 (0.0504)
# M00054 tryptophan
# M01712 cortisol
# M01769 cortisone
# M02342 serotonin (5HT)
# M15140 kynurenine

cd $path_heritability
file="$path_metabolites/metabolites_meta/M02342.metal.pos.txt.gz"
# Extract and organize information from summary.
# Write information to new, temporary file.
echo "identifier allele_1 allele_2 count effect p_value" > summary.txt
zcat $file | awk 'NR > 1 {print $1, $2, $3, $16, $8, $10}' >> summary.txt
head -10 summary.txt
#$path_ldsc/munge_sumstats.py \
#--sumstats summary.txt \
#--out test \
#--merge-alleles $path_alleles/w_hm3.snplist

#$path_ldsc/ldsc.py \
#--h2 $path_access/test.sumstats.gz \
#--ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
#--w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
#--out test_h2

# less test_h2.log
