#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################


# Organize variables.
file_name=$1 # name of file with GWAS summary statistics
path_file=$2 # complete path to file with GWAS summary statistics
path_temporary_gwas_format=$3 # complete path to file for new format
path_temporary_gwas_format_zip=$4 # complete path to file for new format after compression
path_calculate_z_score=$5 # complete path to script to use for z-score standardization
report=$6 # whether to print reports

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "file name: " $file_name
  echo "path to original file: " $path_file
  echo "path to new file: " $path_temporary_gwas_format
  echo "----------"
  echo "----------"
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
rm $path_temporary_gwas_format_zip

# Organize information from linear GWAS.
echo "SNP A1 A2 N BETA P" > $path_temporary_gwas_format
zcat $path_file | awk 'BEGIN { FS=","; OFS=" " } NR > 1 {print $1, toupper($4), toupper($5), $14, $10, $12}' >> $path_temporary_gwas_format
# Calculate Z-score standardization of Beta coefficients.
/usr/bin/bash $path_calculate_z_score \
5 \
$path_temporary_gwas_format \
$path_temporary_gwas_format \
$report

gzip $path_temporary_gwas_format

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "----------"
  echo "----------"
  echo "file name: " $file_name
  head -10 $path_temporary_gwas_format
  echo "----------"
  echo "----------"
  echo "----------"
fi
