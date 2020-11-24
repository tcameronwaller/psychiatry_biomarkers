#!/bin/bash

#chmod u+x script.sh

# Read private, local file paths.
echo "read private file path variables..."
cd ~/paths
path_temporary=$(<"./temporary_bipolar_metabolism.txt")
path_waller="$path_temporary/waller"
path_access="$path_waller/dock/access"
path_disequilibrium="$path_access/disequilibrium"
path_baseline="$path_access/baseline"
path_weights="$path_access/weights"
path_frequencies="$path_access/frequencies"
path_alleles="$path_access/alleles"
path_metabolites="$path_access/metabolites"
path_metabolite_summaries=$(<"./24816252_shin_2014.txt")

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

# Remove previous version of program.
echo "remove previous versions of data..."
rm -r $path_access

# Determine whether the temporary directory structure already exists.
if [ ! -d $path_access ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_access
    mkdir -p $path_disequilibrium
    mkdir -p $path_baseline
    mkdir -p $path_weights
    mkdir -p $path_frequencies
    mkdir -p $path_alleles
    mkdir -p $path_metabolites
fi

cd $path_access
# Linkage disequilibrium scores for European population.
# For simple heritability estimation.
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
tar -xjvf eur_w_ld_chr.tar.bz2 -C $path_disequilibrium
# dock/access/disequilibrium/eur_w_ld_chr/*

# Baseline model linkage disequilibrium scores.
# For partitioned heritability estimation by stratified LD score regression.
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz
tar -xzvf 1000G_Phase3_baselineLD_v2.2_ldscores.tgz -C $path_baseline
# dock/access/baseline/_____/baselineLD.*

# Definitions of Simple Nucleotide Variant alleles.
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
bunzip2 "$path_access/w_hm3.snplist.bz2"
mv "$path_access/w_hm3.snplist" "$path_alleles/w_hm3.snplist"
# w_hm3.snplist

# Weights.
# For partitioned heritability estimation by stratified LD score regression.
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz
tar -xzvf 1000G_Phase3_weights_hm3_no_MHC.tgz -C $path_weights
# dock/access/weights/___/weights.hm3_noMHC.*

# Frequencies.
# For partitioned heritability estimation by stratified LD score regression.
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq.tgz
tar -xzvf 1000G_Phase3_frq.tgz -C $path_frequencies
# dock/access/frequencies/1000G_Phase3_frq/1000G.EUR.QC.*

# GWAS summary statistics for metabolites.
cp -r $path_metabolite_summaries "$path_access/shin_et_al.metal.out.tar.gz"
tar -zxvf "$path_access/shin_et_al.metal.out.tar.gz" -C $path_metabolites
