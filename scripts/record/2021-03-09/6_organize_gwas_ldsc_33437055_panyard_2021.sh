#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize variables.
metabolite=$1 # unique identifier of the metabolite
file_name=$2 # name of file with GWAS summary statistics
path_file=$3 # complete path to file with GWAS summary statistics
path_temporary_gwas_format=$4 # complete path to file for new format
path_destination_parent=$5 # parent directory for files
path_promiscuity_scripts=$6 # complete path to script to use for z-score standardization
report=$7 # whether to print reports

#path_calculate_z_score="$path_promiscuity_scripts/calculate_z_score_column_4_of_5.sh"
path_calculate_z_score="$path_promiscuity_scripts/calculate_z_score_column_5_of_6.sh"

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "metabolite: " $metabolite
  echo "file name: " $file_name
  echo "path to original file: " $path_file
  echo "path to new file: " $path_temporary_gwas_format
  echo "----------"
fi

# Format of GWAS summary statistics for LDSC.
# https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics
# description: ............................ LDSC column
# variant identifier: ........................ "SNP"
# alternate allele (effect allele): .......... "A1"
# reference allele (non-effect allele): ...... "A2"
# sample size: ............................... "N"
# effect (coefficient or odds ratio): ........ "BETA" or "OR"
# probability (p-value): ..................... "P"

# Format of GWAS summary statistics for PRS-CS.
# https://github.com/getian107/PRScs
# description: ............................ PRS-CS column
# variant identifier: ........................ "SNP"
# alternate allele (effect allele): .......... "A1"
# reference allele (non-effect allele): ...... "A2"
# effect (coefficient or odds ratio): ........ "BETA" or "OR"
# probability (p-value): ..................... "P"

# Remove any previous versions of temporary files.
rm $path_temporary_gwas_format

# Define temporary file to collect table's rows.
# Even temporary files need to have names specific to each metabolite.
# During parallel processing, multiple temporary files will exist
# simultaneously.
path_temporary_collection="${path_destination_parent}/temporary_collection_${metabolite}.txt"

# Organize information from linear GWAS.
echo "SNP A1 A2 N BETA P" > $path_temporary_collection
zcat $path_file | awk 'BEGIN { FS=","; OFS=" " } NR > 1 {print $1, toupper($4), toupper($5), $14, $10, $12}' >> $path_temporary_collection
# Calculate Z-score standardization of Beta coefficients.
/usr/bin/bash $path_calculate_z_score \
5 \
$path_temporary_collection \
$path_temporary_gwas_format \
$report

# Compress file format.
# No need in this situation, since each iteration replaces the previous file.
#gzip -cvf $path_temporary_gwas_format > $path_temporary_gwas_format_zip

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "file name: " $file_name
  echo "before standardization:"
  head -10 $path_temporary_collection
  echo "after standardization:"
  head -10 $path_temporary_gwas_format
  echo "----------"
fi

# Remove previous versions of temporary files.
#rm $path_temporary_collection
