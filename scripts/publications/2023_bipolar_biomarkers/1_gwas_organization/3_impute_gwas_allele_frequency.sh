#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 22 February 2023
# Date, last execution: 24 May 2023
################################################################################
# Note

# The purpose of this script is to impute missing frequencies of effect alleles
# in GWAS summary statistics from the "European" superpopulation of the
# 1000 Genomes Project.


################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_source="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/gwas_vcf_process"
path_directory_product="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/gwas_extra_process"
# Files.
path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2023-02-24.tsv"
# Scripts.
path_directory_promiscuity_scripts="${path_directory_process}/promiscuity/scripts"
path_script_impute_gwas_allele_frequency="${path_directory_promiscuity_scripts}/gwas_clean/impute_gwas_allele_frequency.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product

###########################################################################
# Organize parameters.

report="true"

################################################################################
# Execute procedure.

# Most sets of GWAS summary statistics do not need extra processing.
# Copy the GWAS summary statistics from the GWAS2VCF procedure.

cp $path_directory_source/*.txt.gz $path_directory_product

# Impute missing frequencies of effect alleles in GWAS summary statistics.

##########
# 36093044_mathieu_2022

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/36093044_mathieu_2022_hypothyroidism.txt.gz" \
"${path_directory_product}/36093044_mathieu_2022_hypothyroidism.txt.gz" \
$report

##########
# 35459240_said_2022

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/35459240_said_2022_c_reactive_protein.txt.gz" \
"${path_directory_product}/35459240_said_2022_c_reactive_protein.txt.gz" \
$report

##########
# 34017140_mbatchou_2021

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/34017140_mbatchou_2021_albumin.txt.gz" \
"${path_directory_product}/34017140_mbatchou_2021_albumin.txt.gz" \
$report

##########
# 29875488_sun_2018

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/29875488_sun_2018_complement_c4.txt.gz" \
"${path_directory_product}/29875488_sun_2018_complement_c4.txt.gz" \
$report

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/29875488_sun_2018_follitropin_fsh.txt.gz" \
"${path_directory_product}/29875488_sun_2018_follitropin_fsh.txt.gz" \
$report

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/29875488_sun_2018_follistatin.txt.gz" \
"${path_directory_product}/29875488_sun_2018_follistatin.txt.gz" \
$report

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/29875488_sun_2018_lutropin.txt.gz" \
"${path_directory_product}/29875488_sun_2018_lutropin.txt.gz" \
$report

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/29875488_sun_2018_lutropin_beta.txt.gz" \
"${path_directory_product}/29875488_sun_2018_lutropin_beta.txt.gz" \
$report

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/29875488_sun_2018_parathyrin.txt.gz" \
"${path_directory_product}/29875488_sun_2018_parathyrin.txt.gz" \
$report

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/29875488_sun_2018_thyroid_peroxidase.txt.gz" \
"${path_directory_product}/29875488_sun_2018_thyroid_peroxidase.txt.gz" \
$report

##########
# 24586183_medici_2014

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/24586183_medici_2014_thyroid_peroxidase_antibody.txt.gz" \
"${path_directory_product}/24586183_medici_2014_thyroid_peroxidase_antibody.txt.gz" \
$report

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/24586183_medici_2014_thyroid_peroxidase_reactivity.txt.gz" \
"${path_directory_product}/24586183_medici_2014_thyroid_peroxidase_reactivity.txt.gz" \
$report



#
