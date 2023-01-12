#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Review: TCW; 23 November 2022
# Review: TCW; 16 November 2022
# TCW confirmed that the file paths and file names for GWAS summary statistics
# match those in file "table_gwas_summary_heritability_ldsc.ods".

###########################################################################
###########################################################################
###########################################################################



################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")
#path_gwas_summaries_waller=$(<"./gwas_summaries_waller_metabolism.txt")
path_directory_gwas_summaries_team=$(<"./gwas_summaries_team.txt")

path_directory_dock="$path_process/dock"
path_directory_product="${path_directory_dock}/epilepsy/gwas_format_standard"

# Scripts.

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product

# Report.
report="true"

###########################################################################
# Execute procedure.

# Organize multi-dimensional array of information about studies.
studies=()
# [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]

##########

studies+=(
  "${path_directory_gwas_summaries_team}/REFORMATTED_CURRENT/BD_PGC3_EUR.txt.gz;\
  ${path_directory_product}/BD_PGC3_EUR.txt.gz"
)
studies+=(
  "${path_directory_gwas_summaries_team}/REFORMATTED_CURRENT/BDI_PGC3_EUR.txt.gz;\
  ${path_directory_product}/BDI_PGC3_EUR.txt.gz"
)
studies+=(
  "${path_directory_gwas_summaries_team}/REFORMATTED_CURRENT/BDII_PGC3_EUR.txt.gz;\
  ${path_directory_product}/BDII_PGC3_EUR.txt.gz"
)
studies+=(
  "${path_directory_gwas_summaries_team}/REFORMATTED_CURRENT/EPILEPSYgen_ILAEC_TE.txt.gz;\
  ${path_directory_product}/EPILEPSYgen_ILAEC_TE.txt.gz"
)

################################################################################

# Organize information for batch instances.
for study_details in "${studies[@]}"; do
  # Separate fields from instance.
  # [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]
  IFS=";" read -r -a array <<< "${study_details}"
  path_file_source="${array[0]}"
  path_file_product="${array[1]}"
  # Copy the source file path to the product file path.
  cp $path_file_source $path_file_product
  # Report.
  if [[ "$report" == "true" ]]; then
    echo "----------"
    echo "Source file path:"
    echo $path_file_source
    echo "Product file path:"
    echo $path_file_product
    echo "----------"
  fi
done

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "1_access_copy_gwas.sh"
  echo "----------"
fi



#
