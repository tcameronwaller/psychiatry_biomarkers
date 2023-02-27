#!/bin/bash

################################################################################
################################################################################
################################################################################
# Author: T. Cameron Waller
# Date, first execution: 22 February 2023
# Date, last execution: 22 February 2023
################################################################################
################################################################################
################################################################################
# Note

# The purpose of this script is to perform any necessary procedures on GWAS
# summary statistics after completion of the GWAS2VCF procedure.
# 1. Imputation of frequencies of effect alleles


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
#path_directory_product="${path_directory_dock}/hormone_genetics/gwas_vcf_hormones"
# Files.
path_file_translation="${path_directory_parameters}/table_gwas_translation_tcw_2023-02-24.tsv"
# Scripts.
path_promiscuity_scripts="${path_directory_process}/promiscuity/scripts"
path_script_impute_gwas_allele_frequency="${path_promiscuity_scripts}/gwas_clean/impute_gwas_allele_frequency.sh"

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

cp -r $path_directory_source $path_directory_product

# Perform extra procedures on the sets of GWAS summary statistics for which they
# are necessary.

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/36093044_mathieu_2022_hypothyroidism.txt.gz" \
"${path_directory_product}/36093044_mathieu_2022_hypothyroidism.txt.gz" \
$report



#
