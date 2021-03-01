#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# ...
###########################################################################
###########################################################################
###########################################################################

echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "Access summary statistics from GWAS on human metabolome."
echo "Organize files in format suitable for LDSC."
echo "version 1"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo "----------------------------------------------------------------------"
echo ""
echo ""
echo ""

# Organize paths.
# Read private, local file paths.
echo "read private file path variables and organize paths..."
cd ~/paths
path_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_24816252_shin_2014="$path_gwas_summaries/24816252_shin_2014"
path_31959995_schlosser_2021="$path_gwas_summaries/31959995_schlosser_2021"
path_33437055_panyard_2021="$path_gwas_summaries/33437055_panyard_2021"

path_temporary=$(<"./processing_bipolar_metabolism.txt")
path_waller="$path_temporary/waller"
path_bipolar_metabolism="$path_waller/bipolar_metabolism"
path_scripts="$path_waller/bipolar_metabolism/scripts/record/2021-03-02"
path_calculate_z_score="$path_scripts/calculate_z_score.sh"

path_dock="$path_waller/dock"
path_access_gwas="$path_dock/access/gwas"
path_access_shin_2014="$path_access_gwas/24816252_shin_2014"
path_access_schlosser_2021="$path_access_gwas/31959995_schlosser_2021"
path_access_panyard_2021="$path_access_gwas/33437055_panyard_2021"

# Initialize directories.
#rm -r $path_gwas
if [ ! -d $path_access_gwas ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_access_gwas
    mkdir -p $path_access_shin_2014
    mkdir -p $path_access_schlosser_2021
    mkdir -p $path_access_panyard_2021
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
# description: ............................ LDSC column
# variant identifier: ........................ "SNP"
# alternate allele (effect allele): .......... "A1"
# reference allele (non-effect allele): ...... "A2"
# effect (coefficient or odds ratio): ........ "BETA" or "OR"
# probability (p-value): ..................... "P"


# PubMed: 33437055; Author: Panyard; Year: 2021.
echo "----------------------------------------------------------------------"
echo "PubMed: 33437055; Author: Panyard; Year: 2021"
echo "Human genome version: GRCh37, hg19"
echo "----------------------------------------------------------------------"
cd $path_33437055_panyard_2021
metabolite_files=(metabolite_*_meta_analysis_gwas.csv.gz)
count=${#metabolite_files[@]}
echo "count of file paths: " $count
#rm $path_metabolites/metabolite_files.txt
#for path_file in "${metabolite_files[@]}"; do
#    echo $path_file >> $path_metabolites/metabolite_files.txt
#done

# Define glob pattern to recognize relevant files.
pattern="${path_33437055_panyard_2021}/metabolite_*_meta_analysis_gwas.csv.gz"
# Iterate on all files and directories in parent directory.
for file in $path_33437055_panyard_2021/*; do
  if [ -f "$file" ]; then
    # Current content item is a file.
    echo $file
    if [[ "$file" == ${pattern} ]]; then
      # File name matches glob pattern.
      echo "... pattern match! ..."
      base_name="$(basename -- $file)"
      echo "file: " $base_name
      # Copy the file to new directory.
      #cp $file "$path_access_metabolites/$base_name"
    fi
  fi
done





if false; then
  # Organize information from linear GWAS.

  echo "SNP A1 A2 N BETA P" > $path_gwas_alcohol_format
  # Format of GWAS summary for LDSC.
  # https://github.com/bulik/ldsc/wiki/Heritability-and-Genetic-Correlation#reformatting-summary-statistics
  # description: ............................ PGC column ... LDSC column
  # variant identifier: ....................... "SNP" ........ "SNP"
  # alternate allele (effect allele): ......... "A1" ......... "A1"
  # reference allele (non-effect allele): ..... "A2" ......... "A2"
  # sample size: .............................. 52,848 ....... "N"
  # effect (beta): ............................ "Z" .......... "BETA"
  # probability (p-value): .................... "P" .......... "P"

  # SNP: split($2,a,":"); print a[1]
  # A1: toupper($4)
  # A2: toupper($5)
  # N: (cases: 14,904; controls: 37,944; total: 52,848)
  # BETA: $6
  # P: $7
  zcat $path_gwas_alcohol_raw | awk 'NR > 1 {split($2,a,":"); print a[1], toupper($4), toupper($5), (52848), $6, $7}' >> $path_gwas_alcohol_format
  # Calculate Z-score standardization of Beta coefficients.
  #/usr/bin/bash $path_calculate_z_score 5 $path_gwas_alcohol_format $path_gwas_alcohol_format
  gzip $path_gwas_alcohol_format
  echo "after format..."
  head -30 "$path_gwas_alcohol_format.gz"

  echo "----------"
  echo "----------"
  echo "----------"
fi
