#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 2 August 2023
# Date, last execution: 2 August 2023
# Date, review: _ August 2023
################################################################################
# Note

# The purpose of this script is to impute missing frequencies of effect alleles
# in GWAS summary statistics from the "European" superpopulation of the
# 1000 Genomes Project.

# After running this script, check the product directory to make sure that this
# procedure wrote the appropriate files again, later than all others.



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_source="${path_directory_dock}/gwas_disorders_tcw_2023-08-02/2_gwas_assembly_grch37"
path_directory_product="${path_directory_dock}/gwas_disorders_tcw_2023-08-02/3_gwas_allele_frequency"
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
# 30482948_walters_2018

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/30482948_walters_2018_eur_all.txt.gz" \
"${path_directory_product}/30482948_walters_2018_eur_all.txt.gz" \
$report

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/30482948_walters_2018_eur_unrel_meta.txt.gz" \
"${path_directory_product}/30482948_walters_2018_eur_unrel_meta.txt.gz" \
$report

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/30482948_walters_2018_eur_unrel_genotype.txt.gz" \
"${path_directory_product}/30482948_walters_2018_eur_unrel_genotype.txt.gz" \
$report

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/30482948_walters_2018_female.txt.gz" \
"${path_directory_product}/30482948_walters_2018_female.txt.gz" \
$report

/usr/bin/bash $path_script_impute_gwas_allele_frequency \
"${path_directory_source}/30482948_walters_2018_male.txt.gz" \
"${path_directory_product}/30482948_walters_2018_male.txt.gz" \
$report


#
