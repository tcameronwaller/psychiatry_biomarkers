#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 21 September 2023
# Date, last execution: 29 September 2023
# Date, review: 29 September 2023
################################################################################
# Note

# Perform this procedure judiciously.
# The purpose of this procedure is to run a few final checks to prepare GWAS
# summary statistics for analysis in LDSC, LDpred2, or SBayesR.
# 1. Remove any records (rows) with empty values (cells).
# 2. Constrain probabilities (p-values) to double-float precision
# (1E-308 to __).
# 3. Remove records for any SNPs with missing counts of observations, as these
# SNPs raise an error in LDpred2.


# Zhou will need a filter to remove records (rows) with missing count of observations (sample size)...
# either before or after GWAS2VCF.

##########
# TCW; 29 September 2023
# None of the studies below have any missingness or nonsense characters (any
# characters other than "T", "C", "G", or "A") in their designations of effect
# and other alleles.
# 32581359_saevarsdottir_2020_thyroid_autoimmunity
# 32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute
# 32769997_zhou_2020_thyroid_hormone

##########
# TCW; 29 September 2023
# Review of other studies that failed LDpred2 polygenic score procedure.
# These studies failed Anthony Batzler's implementation of LDpred2 polygenic
# score procedure on 26 August 2023.
# 1. Studies that failed for probable reason of weak estimates of SNP
# heritability.
# 29875488_sun_2018_follistatin
# 30367059_teumer_2018_hypothyroidism
# 31169883_pott_2019_aldosterone_all
# 31169883_pott_2019_androstenedione_female
# 31169883_pott_2019_cortisol_female
# 31169883_pott_2019_testosterone_female
# 33441150_dennis_2021_lutropin
# 33441150_dennis_2021_thyroxine_total
# 34822396_pott_2021_aldosterone_all
# 34822396_pott_2021_aldosterone_male
# 34822396_pott_2021_progesterone_all
# 34822396_pott_2021_progesterone_female
# 34822396_pott_2021_testosterone_estradiol_female
# 2. Studies that failed for miscellaneous reasons.
# 36376304_koskeridis_2022_c_reactive_protein        (values of probability greater than zero but less than 1E-307)
# 34662886_backman_2021_albumin                      (low count of total SNPs < 500, 000; low proportion map to HapMap2)
# 34226706_barton_2021_albumin                       (missingness or zero in values of effect and standard error)
# 32769997_zhou_2020_thyroid_hormone                 (missingness in values of count of observations)
# 32581359_saevarsdottir_2020_thyroid_autoimmunity   (unknown reason even after thorough checks)
#   This study originally failed LDpred2 procedure due to genomic coordinates in
#   assembly GRCh38 instead of GRCh37.
#   After correcting the assembly of genomic coordinates, there still seem to be
#   problems with this study that cause errors in GWAS2VCF.
#   https://github.com/pysam-developers/pysam/blob/master/pysam/libcbcf.pyx
#   line 3249
#    if b'' in values:
#     raise ValueError('cannot set null allele')
#   There might be a situation of a reference allele without any alternate alleles,
#   and this situation might arise in the process of filters within GWAS2VCF.
# 32242144_revez_2020_vitamin_d                      (original, raw GWAS summary statistics lacked genomic coordinates for chromosome and base position)

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_parent_temporary="${path_directory_process}/temporary_check_4"
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_source="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-29/3_gwas_allele_frequency"
path_directory_product="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-29/4_filter_constrain_gwas_values"
# Files.
# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_script_check="${path_directory_partner_scripts}/gwas_clean/check_gwas_summary_values.sh"
path_script_filter="${path_directory_partner_scripts}/gwas_clean/filter_constrain_gwas_summary_values.sh"

# Initialize directories.
rm -r $path_directory_parent_temporary # caution
mkdir -p $path_directory_parent_temporary
rm -r $path_directory_product # caution
mkdir -p $path_directory_product
cd $path_directory_product

################################################################################
# Organize parameters.

report="true"

################################################################################
# Execute procedure.

##########
# Copy the GWAS summary statistics from the previous process.
# Most sets of GWAS summary statistics do not need extra processing.
# Subsequent processes on a few studies will replace the appropriate files.
cp $path_directory_source/*.txt.gz $path_directory_product


##########
# 36376304_koskeridis_2022_c_reactive_protein

# records that do not pass checks: 3 of 6,206,409 (header line); most due to values out of range in probability (p-value)
/usr/bin/bash $path_script_check \
"${path_directory_source}/36376304_koskeridis_2022_c_reactive_protein.txt.gz" \
$path_directory_parent_temporary \
$report

if true; then
  # source lines: 6,206,409
  # product lines: 6,206,409
  /usr/bin/bash $path_script_filter \
  "${path_directory_source}/36376304_koskeridis_2022_c_reactive_protein.txt.gz" \
  "${path_directory_product}/36376304_koskeridis_2022_c_reactive_protein.txt.gz" \
  $report

  # records that do not pass checks: 1 of 6,206,409 (header line)
  /usr/bin/bash $path_script_check \
  "${path_directory_product}/36376304_koskeridis_2022_c_reactive_protein.txt.gz" \
  $path_directory_parent_temporary \
  $report
fi



##########
# 34662886_backman_2021_albumin

# records that do not pass checks: 1 of 499,614 (header line)
/usr/bin/bash $path_script_check \
"${path_directory_source}/34662886_backman_2021_albumin.txt.gz" \
$path_directory_parent_temporary \
$report

if true; then
  # source lines: 499,614
  # product lines: 499,614
  /usr/bin/bash $path_script_filter \
  "${path_directory_source}/34662886_backman_2021_albumin.txt.gz" \
  "${path_directory_product}/34662886_backman_2021_albumin.txt.gz" \
  $report
fi



##########
# 34226706_barton_2021_albumin

# records that do not pass checks: 246,495 of 5,515,076 (header line); most due to missingness in effect and standard error
/usr/bin/bash $path_script_check \
"${path_directory_source}/34226706_barton_2021_albumin.txt.gz" \
$path_directory_parent_temporary \
$report

if true; then
  # source lines: 5,515,076
  # product lines: 5,268,582
  /usr/bin/bash $path_script_filter \
  "${path_directory_source}/34226706_barton_2021_albumin.txt.gz" \
  "${path_directory_product}/34226706_barton_2021_albumin.txt.gz" \
  $report
fi



##########
# 32769997_zhou_2020

# records that do not pass checks: 2,993 of 22,397,081 (header line); all due to missingness in count of observations
/usr/bin/bash $path_script_check \
"${path_directory_source}/32769997_zhou_2020_thyroid_hormone.txt.gz" \
$path_directory_parent_temporary \
$report

if true; then
  # source lines: 22,397,081
  # product lines: 22,394,089
  /usr/bin/bash $path_script_filter \
  "${path_directory_source}/32769997_zhou_2020_thyroid_hormone.txt.gz" \
  "${path_directory_product}/32769997_zhou_2020_thyroid_hormone.txt.gz" \
  $report
fi



##########
# 32581359_saevarsdottir_2020

# records that raise checks: 1 of 43,529,207 (header line)
/usr/bin/bash $path_script_check \
"${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz" \
$path_directory_parent_temporary \
$report

if true; then
  # source lines: 43,529,207
  # product lines: 43,529,207
  /usr/bin/bash $path_script_filter \
  "${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz" \
  "${path_directory_product}/32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz" \
  $report
fi

# records that raise checks: 1 of 17,412,123 (header line)
# records with allele frequency of zero: 7,106,273 (header line)
/usr/bin/bash $path_script_check \
"${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz" \
$path_directory_parent_temporary \
$report

if true; then
  # source lines: 17,412,123
  # product lines: 17,412,123
  /usr/bin/bash $path_script_filter \
  "${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz" \
  "${path_directory_product}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz" \
  $report
fi



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "4_check_filter_constrain_gwas_values.sh"
  echo "----------"
fi



#
