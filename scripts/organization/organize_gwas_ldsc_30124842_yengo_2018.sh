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
path_source_file=$3 # full path to source file with GWAS summary statistics
path_genetic_reference=$4 # full path to parent directory with genetic reference files for LDSC
path_gwas=$5 # full path to parent directory for formatted GWAS summary statistics
path_heritability=$6 # full path to parent directory for LDSC heritability estimation
path_genetic_correlation=$7 # full path to parent directory for LDSC genetic correlation
path_promiscuity_scripts=$8 # complete path to directory of scripts for z-score standardization
report=$9 # whether to print reports

################################################################################
# Organize variables.

# Read private, local file paths.
cd ~/paths
path_ldsc=$(<"./tools_ldsc.txt")

path_alleles="$path_genetic_reference/alleles"
path_disequilibrium="$path_genetic_reference/disequilibrium"
path_baseline="$path_genetic_reference/baseline"
path_weights="$path_genetic_reference/weights"
path_frequencies="$path_genetic_reference/frequencies"

path_study_gwas="${path_gwas}/${study}"
path_study_heritability="${path_heritability}/${study}"
path_study_genetic_correlation="${path_genetic_correlation}/${study}"
mkdir -p $path_study_gwas
mkdir -p $path_study_heritability
mkdir -p $path_study_genetic_correlation
path_temporary_collection="${path_study_gwas}/temporary_gwas_collection.txt.gz"
path_gwas_format="${path_study_gwas}/gwas_format.txt.gz"
path_gwas_munge="${path_study_gwas}/gwas_munge"
path_gwas_munge_suffix="${path_gwas_munge}.sumstats.gz"
path_gwas_munge_log="${path_gwas_munge}.log"
path_heritability_report="${path_study_heritability}/heritability_report"

#path_calculate_z_score="$path_promiscuity_scripts/calculate_z_score_column_4_of_5.sh"
path_calculate_z_score="$path_promiscuity_scripts/calculate_z_score_column_5_of_6.sh"

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
  echo "variant identifier (rsID) version: dbSNP151"
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

# Munge GWAS summary statistics for use in LDSC.
$path_ldsc/munge_sumstats.py \
--sumstats $path_gwas_format \
--out $path_gwas_munge \
--merge-alleles $path_alleles/w_hm3.snplist \
#--a1-inc

# Heritability.
$path_ldsc/ldsc.py \
--h2 $path_gwas_munge_suffix \
--ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
--out $path_heritability_report

# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "file name: " $source_file
  echo "before standardization:"
  head -10 $path_temporary_collection
  echo "after standardization:"
  head -10 $path_gwas_format
  echo "after LDSC munge:"
  head -10 $path_gwas_munge_suffix
  echo "LDSC heritability report:"
  cat $path_heritability_report
  echo "----------"
fi

###########################################################################
# Remove temporary files.
rm $path_temporary_collection
