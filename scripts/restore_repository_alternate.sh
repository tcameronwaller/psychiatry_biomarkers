#!/bin/bash

#chmod u+x script.sh

################################################################################
# Note

# The purpose of this script, in contrast to "restore_repository.sh" is to
# pull down alternate versions of the respective repositories without removing
# the current or main versions while they are being used by a running process.

################################################################################
# Organize script parameters.

# Main project that manages and drives execution.
#project="sexy_alcohol"
project="psychiatric_metabolism"

# Read private, local file path for main project directory.
echo "read private file path variables..."
cd ~/paths
path_process=$(<"./process_${project}.txt")
#path_process="/home/tcameronwaller/Downloads"

# Define paths.
path_parameters="$path_process/dock/parameters_alternate" # parent directory for parameter information from multiple repositories
path_repository="$path_process/${project}_alternate" # main project repository
path_partner="$path_process/partner_alternate" # general utility functionality relevant to multiple projects
path_stragglers="$path_process/stragglers_alternate" # functionality specific to smaller studies and data sets
path_uk_biobank="$path_process/uk_biobank_alternate" # functionality specific to data from UK Biobank

# Echo each command to console.
set -x

# Remove previous versions of repositories.

echo "remove previous versions of the repositories..."
rm -r $path_parameters
rm -r $path_repository
rm -r $path_partner
rm -r $path_stragglers
rm -r $path_uk_biobank

##########
# Access and organize current version of the main repository.

echo "access current version of the ${project} repository..."
cd $path_process
wget "https://github.com/tcameronwaller/${project}/archive/main.zip"
unzip main.zip
rm main.zip
mv "${project}-main" $path_repository
mv "$path_repository/package" "$path_repository/${project}"

##########
# Organize and restore supplemental sub-repositories.

# Repository: partner
# Scripts remain within original repository's structure.
# Python code transfers to sub-package.
echo "access current version of the partner repository..."
cd $path_process
wget https://github.com/tcameronwaller/partner/archive/main.zip
unzip main.zip
rm main.zip
mv partner-main $path_partner
mv "$path_partner/package" "$path_partner/partner"
cp -r "$path_partner/partner" "$path_repository/${project}/partner"

# Repository: stragglers
# Scripts remain within original repository's structure.
# Python code transfers to a sub-package within main package.
echo "access current version of the stragglers repository..."
cd $path_process
wget https://github.com/tcameronwaller/stragglers/archive/main.zip
unzip main.zip
rm main.zip
mv stragglers-main $path_stragglers
mv "$path_stragglers/package" "$path_stragglers/stragglers"
cp -r "$path_stragglers/stragglers" "$path_repository/${project}/stragglers"

# Repository: uk_biobank
# Scripts remain within original repository's structure.
# Python code transfers to a sub-package within main package.
echo "access current version of the uk_biobank repository..."
cd $path_process
wget https://github.com/tcameronwaller/uk_biobank/archive/main.zip
unzip main.zip
rm main.zip
mv uk_biobank-main $path_uk_biobank
mv "$path_uk_biobank/package" "$path_uk_biobank/uk_biobank"
cp -r "$path_uk_biobank/uk_biobank" "$path_repository/${project}/uk_biobank"

##########
# Organize and restore parameters.

mkdir -p $path_parameters
cp -r "$path_repository/parameters" "$path_parameters/parameters"
mv "$path_parameters/parameters" "$path_parameters/${project}"
cp -r "$path_partner/parameters" "$path_parameters/parameters"
mv "$path_parameters/parameters" "$path_parameters/partner"
cp -r "$path_stragglers/parameters" "$path_parameters/parameters"
mv "$path_parameters/parameters" "$path_parameters/stragglers"
cp -r "$path_uk_biobank/parameters" "$path_parameters/parameters"
mv "$path_parameters/parameters" "$path_parameters/uk_biobank"
