#!/bin/bash

#chmod u+x script.sh

###########################################################################
###########################################################################
###########################################################################
# This script organizes directories and iteration instances then submits
# script "regress_metabolite_heritability.sh" to the Sun Grid Engine.
###########################################################################
###########################################################################
###########################################################################

echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "----------"
echo "The script organizes and submits an array batch job."
echo "----------"
echo "version check: 1"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"

# Organize paths.
# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths
path_temporary=$(<"./temporary_bipolar_metabolism.txt")
path_waller="$path_temporary/waller"
path_dock="$path_waller/dock"
path_disequilibrium="$path_dock/access/disequilibrium"
path_alleles="$path_dock/access/alleles"
path_metabolites="$path_dock/access/metabolites"
path_metabolite_summaries="$path_dock/access/metabolites/metabolites_meta"
path_heritability="$path_dock/heritability"
path_heritability_metabolites="$path_dock/heritability/metabolites"
path_bipolar_metabolism="$path_waller/bipolar_metabolism"
path_scripts="$path_bipolar_metabolism/scripts/metabolite_heritability"

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

# Determine whether the temporary directory structure already exists.
if [ ! -d $path_heritability_metabolites ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_heritability_metabolites
fi

# Organize instances for iteration.
echo "----------------------------------------------------------------------"
echo "Organize the array of batch instances."
echo "----------------------------------------------------------------------"
cd $path_metabolite_summaries
echo $path_metabolite_summaries
metabolite_files=(*.metal.pos.txt.gz)
#printf "%s\n" "${metabolite_files[@]}" > $path_metabolites/metabolite_files.txt
count=${#metabolite_files[@]}
echo "count of file paths: " $count
#cat $metabolite_files
for path_file in "${metabolite_files[@]}"; do
    echo $ path_file >> $path_metabolites/metabolite_files.txt
done


# Submit array batch to Sun Grid Engine.
# Array batch indices cannot start at zero.
echo "----------------------------------------------------------------------"
echo "Submit array batch to Sun Grid Engine."
echo "----------------------------------------------------------------------"
#qsub -t 1-${count}:1 -o \
#qsub -t 1-5:1 -o \
#"$path_heritability/out.txt" -e "$path_heritability/error.txt" \
#$path_scripts/regress_metabolite_heritability.sh \
#$path_metabolites/metabolite_files.txt \
#$path_dock $count
