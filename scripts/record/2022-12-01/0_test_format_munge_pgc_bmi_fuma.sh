#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Review: TCW; __ November 2022

# Notes; TCW; 18 November 2022:
# 1. LDSC Munge did interpret successfully the SNPs with original identifiers.
# - - "   print $1, $2, $3, toupper($4), toupper($5), "NA", $7, $8, $9, (3717), "NA", (1), "NA", "NA"   "
# 2. Original SNP identifier format was "chr[chromosome]_[base-pair position]_[ambiguous allele]".
# 3. LDSC Munge assigned rsIDs to these SNPs.
# 4. LDSC kept 1,217,311 SNPs (1,144,121 non-missing) from interpretation of original identifiers.
# 5. LDSC Munge did not interpret the SNPs with novel identifiers in format "[chromosome]:[base-pair position]".
# - - "   split($1, a, "_"); (b = a[1]); sub(/chr/, "", b); print (b ":" a[2]), $2, $3, toupper($4), toupper($5), "NA", $7, $8, $9, (3717), "NA", (1), "NA", "NA"   "
# 6. LDSC Munge did not interpret the SNPs with novel identifiers in format "[chromosome]_[base-pair position]".
# - - "   split($1, a, "_"); (b = a[1]); sub(/chr/, "", b); print (b "_" a[2]), $2, $3, toupper($4), toupper($5), "NA", $7, $8, $9, (3717), "NA", (1), "NA", "NA"   "
# 7. LDSC Munge did not interpret the SNPs with novel identifiers in format "[chromosome]_[base-pair position]_[ambiguous allele]".
# - - "   split($1, a, "_"); (b = a[1]); sub(/chr/, "", b); print (b "_" a[2] "_" a[3]), $2, $3, toupper($4), toupper($5), "NA", $7, $8, $9, (3717), "NA", (1), "NA", "NA"   "

###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize paths.

# Read private, local file paths.
cd ~/paths
path_ldsc=$(<"./tools_ldsc.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")

path_dock="$path_process/dock"
path_directory_source="${path_dock}/bipolar_body/gwas_access"
path_directory_product="${path_dock}/bipolar_body/test_format_munge_pgc_bmi"
path_file_gwas_source="${path_directory_source}/bmi_bipolar_case_pgc_mafe_fuma.txt.gz"
path_file_gwas_source_decompress="${path_directory_source}/bmi_bipolar_case_pgc_mafe_fuma.txt"
path_file_gwas_format_team="${path_directory_product}/gwas_format_team.txt"
path_file_gwas_format_ldsc="${path_directory_product}/gwas_format_ldsc.txt"
path_file_base_gwas_munge="${path_directory_product}/gwas_munge_ldsc"
path_directory_reference="${path_dock}/bipolar_body/reference_ldsc"
path_file_alleles="${path_directory_reference}/alleles/w_hm3.snplist"

# Scripts.

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product

###########################################################################
# Execute procedure.

gunzip -cvf $path_file_gwas_source > $path_file_gwas_source_decompress
head $path_file_gwas_source_decompress
wc -l $path_file_gwas_source_decompress

if true; then
  # One-step format directly to LDSC.
  echo "SNP A1 A2 N BETA P" > $path_file_gwas_format_ldsc
  zcat $path_file_gwas_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {print $1, toupper($4), toupper($5), (3717), $7, $9}' >> $path_file_gwas_format_ldsc
  head $path_file_gwas_format_ldsc
else
  # Translation to Team format.
  echo "SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT" > $path_file_gwas_format_team
  zcat $path_file_gwas_source | awk 'BEGIN {FS = " "; OFS = " "} NR > 1 {
    print $1, $2, $3, toupper($4), toupper($5), "NA", $7, $8, $9, (3717), "NA", (1), "NA", "NA"
  }' >> $path_file_gwas_format_team
  head $path_file_gwas_format_team

  # Translation to LDSC format.
  echo "SNP A1 A2 N BETA P" > $path_file_gwas_format_ldsc
  cat $path_file_gwas_format_team | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {
    print $1, toupper($4), toupper($5), $10, $7, $9
  }' >> $path_file_gwas_format_ldsc
  head $path_file_gwas_format_ldsc
fi

# Munge GWAS summary statistics in LDSC.
$path_ldsc/munge_sumstats.py \
--sumstats $path_file_gwas_format_ldsc \
--signed-sumstats BETA,0 \
--merge-alleles $path_file_alleles \
--out $path_file_base_gwas_munge \


#
