#!/bin/bash

#chmod u+x script.sh

###########################################################################
###########################################################################
###########################################################################
# This script organizes paths and parameters to access local information
# about persons and phenotypes in the UKBiobank.
###########################################################################
###########################################################################
###########################################################################

# version: 1

###########################################################################
# Organize script parameters.
#project="sexy_alcohol"
project="psychiatric_metabolism"

# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths
path_scripts_team_uk_biobank=$(<"./ukbiobank_scripts.txt")
path_process=$(<"./process_${project}.txt")

path_scripts_uk_biobank="$path_process/uk_biobank/scripts"
path_script_batzler_access="$path_scripts_team_uk_biobank/pullUKBclinical.sh"

path_dock="${path_process}/dock"
path_import="${path_dock}/access/ukbiobank_import"
path_variables="${path_import}/variables.txt"

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

# Initialize directories.
rm -r $path_import
mkdir -p $path_import

# Specify UKBiobank phenotype variables of interest.
# Only use this method if not using the file in the parameter directory.
#echo "specify phenotype variables to access from UKBiobank"
echo "22009 22006 31 22001 21022 21002 50 21001 23104 30890" | tr -s " " "\n" > $path_variables

# Access phenotype variables and auxiliary information from UKBiobank.
/usr/bin/bash "$path_script_batzler_access" \
-n "waller_import" \
-f $path_variables \
-o $path_import \
-d "defineUKB_variables.R" \
