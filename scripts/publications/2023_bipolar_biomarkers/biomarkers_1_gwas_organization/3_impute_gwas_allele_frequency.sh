#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 22 February 2023
# Date, last execution: 26 September 2023
# Date, review: 26 September 2023
################################################################################
# Note

# The purpose of this script is to impute missing frequencies of effect alleles
# in GWAS summary statistics from the "European" superpopulation of the
# 1000 Genomes Project.

# After running this script, check the product directory to make sure that this
# procedure wrote the appropriate files again, later than all others.

# Make a copy of "Saevarsdottir 2020" and impute European allele frequencies
# to test in case there were problems in the original allele frequencies.



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_source="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-25/2_gwas_assembly_grch37"
path_directory_product="${path_directory_dock}/gwas_biomarkers_tcw_2023-09-25/3_gwas_allele_frequency"
# Files.
# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_script_impute_gwas_allele_frequency="${path_directory_partner_scripts}/gwas_clean/impute_gwas_allele_frequency.sh"

# Initialize directories.
rm -r $path_directory_product # caution
mkdir -p $path_directory_product
cd $path_directory_product

###########################################################################
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
# Impute missing frequencies of effect alleles in GWAS summary statistics.

##########
# 36093044_mathieu_2022

# source lines: 10,836,151
# product lines: 4,053,261
/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/36093044_mathieu_2022_hypothyroidism.txt.gz" \
"${path_directory_product}/36093044_mathieu_2022_hypothyroidism.txt.gz" \
$report

##########
# 35459240_said_2022

# source lines: 11,106,738
# product lines: 3,644,730
/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/35459240_said_2022_c_reactive_protein.txt.gz" \
"${path_directory_product}/35459240_said_2022_c_reactive_protein.txt.gz" \
$report

##########
# 34017140_mbatchou_2021

# source lines: 11,346,500
# product lines: 6,500,544
/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/34017140_mbatchou_2021_albumin.txt.gz" \
"${path_directory_product}/34017140_mbatchou_2021_albumin.txt.gz" \
$report

##########
# 32581359_saevarsdottir_2020

# source lines: 43,529,207
# product lines: 17,412,123
/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz" \
"${path_directory_product}/32581359_saevarsdottir_2020_thyroid_autoimmunity_af_impute.txt.gz" \
$report

##########
# 29875488_sun_2018

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/29875488_sun_2018_complement_c4.txt.gz" \
"${path_directory_product}/29875488_sun_2018_complement_c4.txt.gz" \
$report

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/29875488_sun_2018_follitropin.txt.gz" \
"${path_directory_product}/29875488_sun_2018_follitropin.txt.gz" \
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

# source lines: 2,424,987
# product lines: 727,787
/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/24586183_medici_2014_thyroid_peroxidase_antibody.txt.gz" \
"${path_directory_product}/24586183_medici_2014_thyroid_peroxidase_antibody.txt.gz" \
$report

# source lines: 2,424,987
# product lines: 727,787
/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/24586183_medici_2014_thyroid_peroxidase_reactivity.txt.gz" \
"${path_directory_product}/24586183_medici_2014_thyroid_peroxidase_reactivity.txt.gz" \
$report



#
