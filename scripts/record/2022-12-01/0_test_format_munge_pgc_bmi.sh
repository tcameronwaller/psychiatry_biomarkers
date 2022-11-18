#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Review: TCW; __ November 2022

# next script in pipeline
# 1. translate format of all GWAS summary stats files to the common team format
# 2. translate format of all GWAS summary stats files to the LDSC format

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
path_file_gwas_format="${path_directory_product}/gwas_format_ldsc.txt"
path_file_base_gwas_munge="${path_directory_product}/gwas_munge_ldsc"
path_directory_reference="${path_dock}/bipolar_body/reference_ldsc"
path_file_alleles="${path_directory_reference}/alleles/w_hm3.snplist"

# Scripts.

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product

###########################################################################
# Execute procedure.

echo "SNP A1 A2 N BETA P" > $path_file_gwas_format
zcat $path_file_gwas_source | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {print $1, toupper($4), toupper($5), 4332, $7, $9}' >> $path_file_gwas_format
head $path_file_gwas_format

#--signed-sumstats BETA,0 \ # I don't think this argument is consequential.
$path_ldsc/munge_sumstats.py \
--sumstats $path_file_gwas_format \
--signed-sumstats BETA,0 \
--merge-alleles $path_file_alleles \
--out $path_file_base_gwas_munge \


#
