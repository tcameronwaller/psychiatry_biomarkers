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
echo "read private file path variables..."
cd ~/paths
path_process=$(<"./process_${project}.txt")
path_repository="$path_process/${project}"
path_uk_biobank="$path_process/uk_biobank"
path_promiscuity="$path_process/promiscuity"
path_parameters="$path_process/dock/parameters"
path_scripts="$path_uk_biobank/scripts"

path_variables="$path_parameters/${project}/uk_biobank_access_variables.txt"
path_dock="$path_process/dock"
path_access="$path_dock/access/ukbiobank_phenotypes"

# Echo each command to console.
#set -x
# Suppress echo each command to console.
set +x

# Determine whether the temporary directory structure already exists.
if [ ! -d $path_access ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_access
fi

# Specify UKBiobank phenotype variables of interest.
# Only use this method if not using the file in the parameter directory.
#echo "specify phenotype variables to access from UKBiobank"
#echo "31 22001 21022 21002 30600" | tr -s " " "\n" > $path_variables

# Access phenotype variables and auxiliary information from UKBiobank.
/usr/bin/bash "$path_scripts/access_ukbiobank_phenotypes.sh" \
$path_variables \
$path_access
