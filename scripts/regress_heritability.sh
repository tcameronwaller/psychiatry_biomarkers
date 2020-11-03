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
# M00054 tryptophan
# M01712 cortisol
# M01769 cortisone
# M02342 serotonin (5HT)
# M15140 kynurenine

cd $path_heritability
#gunzip $path_metabolites/metabolites_meta/M02342.metal.pos.txt.gz -c >> raw
zcat $path_metabolites/metabolites_meta/M02342.metal.pos.txt.gz >> raw
$pathldsc/munge_sumstats.py \
--sumstats raw \
--out test \
--merge-alleles $path_alleles/w_hm3.snplist

$pathldsc/ldsc.py \
--h2 $path_access/test.sumstats.gz \
--ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--out test_h2

# less test_h2.log
