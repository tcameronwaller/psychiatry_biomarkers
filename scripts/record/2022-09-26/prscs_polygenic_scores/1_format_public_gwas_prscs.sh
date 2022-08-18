#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Test procedure or rescue.
###########################################################################
###########################################################################
###########################################################################


################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="$path_process/dock"

path_gwas_source="${path_gwas_summaries}/34002096_mullins_2021/daner_bip_pgc3_nm_nomay1.gz"
path_gwas_product_container="${path_dock}/gwas_format_prscs/34002096_mullins_2021_all_no_mayo"

###########################################################################
# Execute procedure.

##############################################################################
# Format GWAS summary statistics for analysis in PRS-CS.

# Initialize directory.
rm -r $path_gwas_product_container
mkdir -p $path_gwas_product_container

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_format_gwas="${path_promiscuity_scripts}/gwas_process/format_gwas_prscs/format_gwas_prscs_34002096_mullins_2021_all_no_mayo.sh"

################################################################################
# Paths.
path_gwas_constraint="${path_gwas_product_container}/gwas_constraint.txt"
path_gwas_format="${path_gwas_product_container}/gwas_format.txt"
path_gwas_format_compress="${path_gwas_product_container}/gwas_format.txt.gz"

################################################################################
# Format adaptation.
report="true" # "true" or "false"
/usr/bin/bash "$path_script_format_gwas" \
$path_gwas_source \
$path_gwas_constraint \
$path_gwas_format \
$path_gwas_format_compress \
$report

#
