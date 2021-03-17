#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize arguments.
study=$1 # identifier of GWAS study
source_file=$2 # name of source file with GWAS summary statistics
path_source_file=$3 # complete path to source file with GWAS summary statistics
path_gwas_format=$4 # full path to file for new format
path_product_directory=$5 # full path to parent directory for product files
path_promiscuity_scripts=$6 # complete path to directory of scripts for z-score standardization
report=$7 # whether to print reports

################################################################################
# Organize variables.

#path_calculate_z_score="$path_promiscuity_scripts/calculate_z_score_column_4_of_5.sh"
path_calculate_z_score="$path_promiscuity_scripts/calculate_z_score_column_5_of_6.sh"
path_temporary_collection="${path_product_directory}/temporary_gwas_collection.txt.gz"

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------------------------------------------------------------------"
  echo "----------------------------------------------------------------------"
  echo "----------------------------------------------------------------------"
  echo "Organize GWAS summary statistics."
  echo "PubMed: 30124842"
  echo "author: Yengo"
  echo "date: 16 August 2018"
  echo "phenotype: body mass index"
  echo "Human genome version: GRCh37, hg19"
  echo "----------------------------------------------------------------------"
  echo "----------------------------------------------------------------------"
  echo "----------------------------------------------------------------------"
  echo ""
  echo ""
  echo ""
  echo "----------"
  echo "file name: " $source_file
  echo "path to original file: " $path_source_file
  echo "path to new file: " $path_gwas_format
  echo "----------"
fi

# Format of GWAS summary statistics for LDSC.
# https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics
# description: ............................ LDSC column ........... source column ..............
# variant identifier (RS ID): .............  "SNP" ................  "SNP" .....................
# alternate allele (effect allele): .......  "A1" .................  "Tested_Allele" ...........
# reference allele (non-effect allele): ...  "A2" .................  "Other_Allele" ............
# sample size: ............................  "N" ..................  "N" .......................
# effect (coefficient or odds ratio): .....  "BETA" or "OR" .......  "BETA" or "BETA_COJO" .....
# probability (p-value): ..................  "P" ..................  "P" or "P_COJO" ...........

# Remove any previous versions of temporary files.
rm $path_termporary_collection
rm $path_gwas_format

# Organize information from linear GWAS.
echo "SNP A1 A2 N BETA P" > $path_temporary_collection
zcat $path_source_file | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {print $3, toupper($4), toupper($5), $10, $7, $9}' >> $path_temporary_collection
# Calculate Z-score standardization of Beta coefficients.
/usr/bin/bash $path_calculate_z_score \
5 \
$path_temporary_collection \
$path_gwas_format \
$report

# Compress file format.
# No need in this situation, since each iteration replaces the previous file.
#gzip -cvf $path_temporary_gwas_format > $path_temporary_gwas_format_zip

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "file name: " $source_file
  echo "before standardization:"
  head -10 $path_temporary_collection
  echo "after standardization:"
  head -10 $path_gwas_format
  echo "----------"
fi

###########################################################################
# Remove temporary files.
rm $path_temporary_collection
