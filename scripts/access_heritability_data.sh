#!/bin/bash

#chmod u+x script.sh

# Read private, local file paths.
echo "read private file path variables..."
cd ~/paths
path_temporary=$(<"./temporary_bipolar_metabolism.txt")
path_waller="$path_temporary/waller"
path_access="$path_waller/dock/access"
path_disequilibrium="$path_access/disequilibrium"
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
    mkdir -p $path_metabolites
fi

# Linkage disequilibrium scores for European population.
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
tar -jxvf eur_w_ld_chr.tar.bz2 -C $path_access

# GWAS summary statistics for metabolites.
cp -r $path_metabolite_summaries "$path_access/shin_et_al.metal.out.tar.gz"
tar -zxvf "$path_access/shin_et_al.metal.out.tar.gz" -C $path_metabolites
