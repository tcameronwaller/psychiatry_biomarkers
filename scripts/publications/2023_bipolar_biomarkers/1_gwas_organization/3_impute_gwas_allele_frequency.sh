#!/bin/bash

################################################################################
################################################################################
################################################################################
# Author: T. Cameron Waller
# Date, first execution: 22 February 2023
# Date, last execution: 27 February 2023
################################################################################
################################################################################
################################################################################
# Note

# The purpose of this script is to perform any necessary procedures on GWAS
# summary statistics after completion of the GWAS2VCF procedure.
# 1. Imputation of frequencies of effect alleles

# Need imputation of allele frequencies
# 36093044_mathieu_2022
# 35459240_said_2022
# 34017140_mbatchou_2021
# 29875488_sun_2018
# 24586183_medici_2014

# Do not need imputation of allele frequencies (was just testing)
# 30367059_teumer_2018


################################################################################
################################################################################
################################################################################



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

###########################################################################
# Execute procedure.
###########################################################################

# Most sets of GWAS summary statistics do not need extra processing.
# Copy the GWAS summary statistics from the GWAS2VCF procedure.

cp $path_directory_source/*.txt.gz $path_directory_product

# Perform extra procedures on the sets of GWAS summary statistics for which they
# are necessary.

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/36093044_mathieu_2022_hypothyroidism.txt.gz" \
"${path_directory_product}/36093044_mathieu_2022_hypothyroidism.txt.gz" \
$report

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/30367059_teumer_2018_hypothyroidism.txt.gz" \
"${path_directory_product}/30367059_teumer_2018_hypothyroidism.txt.gz" \
$report

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/30367059_teumer_2018_hyperthyroidism.txt.gz" \
"${path_directory_product}/30367059_teumer_2018_hyperthyroidism.txt.gz" \
$report


#
