#!/bin/bash

#chmod u+x script.sh

# Read private, local file paths.
echo "read private file path variables..."
cd ~/paths
path_temporary=$(<"./temporary_bipolar_metabolism.txt")
path_waller="$path_temporary/waller"
path_bipolar_metabolism="$path_waller/bipolar_metabolism"
path_promiscuity="$path_waller/promiscuity"

# Echo each command to console.
set -x

# Remove previous version of program.

echo "remove previous versions of the repositories..."
rm -r $path_bipolar_metabolism
rm -r $path_promiscuity

# Access current version of the program.

echo "access current version of the bipolar_metabolism repository..."
cd $path_waller
wget https://github.com/tcameronwaller/bipolar_metabolism/archive/main.zip
unzip main.zip
rm main.zip
mv bipolar_metabolism-main $path_bipolar_metabolism
mv "$path_bipolar_metabolism/package" "$path_bipolar_metabolism/bipolar_metabolism"

echo "access current version of the promiscuity repository..."
cd $path_waller
wget https://github.com/tcameronwaller/promiscuity/archive/main.zip
unzip main.zip
rm main.zip
mv promiscuity-main $path_promiscuity
mv "$path_promiscuity/package" "$path_promiscuity/promiscuity"

echo "organize promiscuity as a subpackage within bipolar_metabolism..."
cp -r "$path_promiscuity/promiscuity" "$path_bipolar_metabolism/bipolar_metabolism/promiscuity"
