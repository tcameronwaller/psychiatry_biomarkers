#!/bin/bash

###########################################################################
# Organize paths.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths

path_source_mayo_bipolar_genotype=$(<"./mayo_bipolar_genotype.txt")

path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="${path_process}/dock"
path_target_container="${path_dock}/access"
path_target_mayo_bipolar_genotype="${path_dock}/access/mayo_bipolar_genotype_raw"

###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
set -x

cp -r "$path_source_mayo_bipolar_genotype" "$path_target_container"
mv "$path_target_container/MERGED" "$path_target_mayo_bipolar_genotype"
