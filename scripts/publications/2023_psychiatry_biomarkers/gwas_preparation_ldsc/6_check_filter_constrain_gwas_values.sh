#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 2 October 2023
# Date, last execution: 2 October 2023
# Date, review: 2 October 2023
################################################################################
# Note

# Temporary note about sequence of steps...
# 1_translate_gwas_to_standard_format.sh
# 2_translate_gwas_assembly_to_grch37.sh
# 3_fill_dummy_gwas_allele_frequency.sh
# 4_check_filter_constrain_gwas_values.sh
# 5_clean_gwas_gwas2vcf.sh
# 6_check_filter_constrain_gwas_values.sh

# Study: 32581359_saevarsdottir_2020_thyroid_autoimmunity
# This study failed GWAS2VCF procedure for unknown reasons.
# Rescue this study by copying set of GWAS summary statistics from step 4 of
# preparation procedure before GWAS2VCF.

# Study: 30718901_howard_2019_pgc_ukb
# This study failed GWAS2VCF procedure because the original GWAS summary
# statistics did not include information about chromosome and base-pair
# coordinates. This information is not necessary for analysis in LDSC.
# Rescue this study by copying set of GWAS summary statistics from step 4 of
# preparation procedure before GWAS2VCF.

# Study: 32242144_revez_2020_vitamin_d
# This study failed GWAS2VCF procedure due to missing information in the raw
# GWAS summary statistics.
# Remove this study from the collection.
# It will be necessary to accommodate this exclusion in subsequent procedures.

################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_parent_temporary="${path_directory_process}/temporary_check_6"
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"

path_directory_source="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-29/5_gwas_clean_gwas2vcf"
path_directory_source_rescue="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-29/4_filter_constrain_gwas_values"
path_directory_product="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-29/6_filter_constrain_gwas_values"

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
# 36376304_koskeridis_2022_c_reactive_protein

# records that do not pass checks: 1 of 6,205,745 (header line)
/usr/bin/bash $path_script_check \
"${path_directory_source}/36376304_koskeridis_2022_c_reactive_protein.txt.gz" \
$path_directory_parent_temporary \
$report

##########
# 34662886_backman_2021_albumin

# records that do not pass checks: 1 of 444,842 (header line)
/usr/bin/bash $path_script_check \
"${path_directory_source}/34662886_backman_2021_albumin.txt.gz" \
$path_directory_parent_temporary \
$report

##########
# 34226706_barton_2021_albumin

# records that do not pass checks: 1 of 4,746,640 (header line)
/usr/bin/bash $path_script_check \
"${path_directory_source}/34226706_barton_2021_albumin.txt.gz" \
$path_directory_parent_temporary \
$report

##########
# 32769997_zhou_2020

# records that do not pass checks: 1 of 21,795,381 (header line)
/usr/bin/bash $path_script_check \
"${path_directory_source}/32769997_zhou_2020_thyroid_hormone.txt.gz" \
$path_directory_parent_temporary \
$report

##########
# 32581359_saevarsdottir_2020

# records that raise checks: 1 of 1 (header line)
/usr/bin/bash $path_script_check \
"${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz" \
$path_directory_parent_temporary \
$report

# records that raise checks: 1 of 17,411,881 (header line)
/usr/bin/bash $path_script_check \
"${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz" \
$path_directory_parent_temporary \
$report



##########
# Copy the GWAS summary statistics from the previous process.
# Most sets of GWAS summary statistics do not need extra processing.
# Subsequent processes on a few studies will replace the appropriate files.
cp $path_directory_source/*.txt.gz $path_directory_product



##########
# Remove files of GWAS summary statistics that failed GWAS2VCF procedure.

rm $path_directory_product/32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz
rm $path_directory_product/30718901_howard_2019_pgc_ukb.txt.gz
rm $path_directory_product/32242144_revez_2020_vitamin_d.txt.gz



##########
# Copy GWAS summary statistics from before GWAS2VCF where appropriate.

cp $path_directory_source_rescue/32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz $path_directory_product
cp $path_directory_source_rescue/30718901_howard_2019_pgc_ukb.txt.gz $path_directory_product



################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "6_check_filter_constrain_gwas_values.sh"
  echo "----------"
fi



#
