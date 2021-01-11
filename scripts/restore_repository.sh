#!/bin/bash

#chmod u+x script.sh

###########################################################################
# Organize script parameters.
#project="sexy_alcohol"
project="bipolar_metabolism"

# Read private, local file paths.
echo "read private file path variables..."
cd ~/paths
path_temporary=$(<"./processing_${project}.txt")
path_repository="$path_temporary/waller/${project}"
path_uk_biobank="$path_temporary/waller/uk_biobank"
path_promiscuity="$path_temporary/waller/promiscuity"
path_parameters="$path_temporary/waller/dock/parameters"

# Echo each command to console.
set -x

# Remove previous version of program.

echo "remove previous versions of the repositories..."
rm -r $path_repository
rm -r $path_uk_biobank
rm -r $path_promiscuity
rm -r $path_parameters

##########
# Access and organize current version of the main repository.

echo "access current version of the ${project} repository..."
cd $path_waller
wget "https://github.com/tcameronwaller/${project}/archive/main.zip"
unzip main.zip
rm main.zip
mv "${project}-main" $path_repository
mv "$path_repository/package" "$path_repository/${project}"

##########
# Organize and restore supplemental sub-repositories.

# Repository: uk_biobank
# Scripts remain within original repository's structure.
# Python code transfers to a sub-package within main package.
echo "access current version of the uk_biobank repository..."
cd $path_waller
wget https://github.com/tcameronwaller/uk_biobank/archive/main.zip
unzip main.zip
rm main.zip
mv uk_biobank-main $path_uk_biobank
mv "$path_uk_biobank/package" "$path_uk_biobank/uk_biobank"
cp -r "$path_uk_biobank/uk_biobank" "$path_repository/${project}/uk_biobank"

# Repository: promiscuity
# Scripts remain within original repository's structure.
# Python code transfers to sub-package.
echo "access current version of the promiscuity repository..."
cd $path_waller
wget https://github.com/tcameronwaller/promiscuity/archive/main.zip
unzip main.zip
rm main.zip
mv promiscuity-main $path_promiscuity
mv "$path_promiscuity/package" "$path_promiscuity/promiscuity"
cp -r "$path_promiscuity/promiscuity" "$path_repository/${project}/promiscuity"

##########
# Organize and restore parameters.

mkdir -p $path_parameters
cp -r "$path_repository/parameters" "$path_parameters/parameters"
mv "$path_parameters/parameters" "$path_parameters/${project}"
cp -r "$path_uk_biobank/parameters" "$path_parameters/parameters"
mv "$path_parameters/parameters" "$path_parameters/uk_biobank"
cp -r "$path_promiscuity/parameters" "$path_parameters/parameters"
mv "$path_parameters/parameters" "$path_parameters/promiscuity"
