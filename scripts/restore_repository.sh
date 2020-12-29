#!/bin/bash

#chmod u+x script.sh

# Read private, local file paths.
echo "read private file path variables..."
cd ~/paths
path_temporary=$(<"./temporary_bipolar_metabolism.txt")
path_waller="$path_temporary/waller"
path_bipolar_metabolism="$path_waller/bipolar_metabolism"
path_uk_biobank="$path_waller/uk_biobank"
path_promiscuity="$path_waller/promiscuity"

# Echo each command to console.
set -x

# Remove previous version of program.

echo "remove previous versions of the repositories..."
rm -r $path_bipolar_metabolism
rm -r $path_uk_biobank
rm -r $path_promiscuity

##########
# Access and organize current version of the main repository.

echo "access current version of the bipolar_metabolism repository..."
cd $path_waller
wget https://github.com/tcameronwaller/bipolar_metabolism/archive/main.zip
unzip main.zip
rm main.zip
mv bipolar_metabolism-main $path_bipolar_metabolism
mv "$path_bipolar_metabolism/package" "$path_bipolar_metabolism/bipolar_metabolism"



##########
# Organize and restore parameters.

# TODO: copy the "parameters" directory from the repository to the "dock"...
# TODO: modify the access procedure to copy this parameters information to dock, also...
# Transfer new version of parameters.
#rm $path_access_table_variables
#cp $path_parameters_table_variables $path_access_table_variables

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
cp -r "$path_uk_biobank/uk_biobank" "$path_bipolar_metabolism/bipolar_metabolism/uk_biobank"

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
cp -r "$path_promiscuity/promiscuity" "$path_bipolar_metabolism/bipolar_metabolism/promiscuity"
