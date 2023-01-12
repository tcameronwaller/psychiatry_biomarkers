#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Review: TCW; __ November 2022

###########################################################################
###########################################################################
###########################################################################

################################################################################
# Organize paths.

# Read private, local file paths.
cd ~/paths
path_directory_process=$(<"./process_psychiatric_metabolism.txt")

path_directory_dock="$path_directory_process/dock"
path_directory_source="${path_directory_dock}/epilepsy/gwas_munge_ldsc"
path_directory_product="${path_directory_dock}/epilepsy/gwas_genetic_correlation_ldsc"
path_directory_reference="${path_directory_dock}/hormone_genetics/reference_ldsc"
path_directory_disequilibrium="${path_directory_reference}/disequilibrium/eur_w_ld_chr"

# Scripts.
path_directory_promiscuity_scripts="${path_directory_process}/promiscuity/scripts"
path_file_script="${path_directory_promiscuity_scripts}/utility/ldsc/estimate_gwas_genetic_correlation_ldsc.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product

################################################################################
# Organize parameters.

##########
# Common parameters.
threads=8
report="true"

##########
# Organize multi-dimensional array of information about studies.
# [full path to base name of product file] ; \
# [full path to primary source file of LDSC munge GWAS summary statistics] ; \
# [full path to secondary source file of LDSC munge GWAS summary statistics]

comparisons=()

comparisons+=(
  "${path_directory_product}/bipolar_pgc3_against_epilepsy;\
  ${path_directory_source}/BD_PGC3_EUR.sumstats.gz;\
  ${path_directory_source}/EPILEPSYgen_ILAEC_TE.sumstats.gz"
)
comparisons+=(
  "${path_directory_product}/bipolar1_pgc3_against_epilepsy;\
  ${path_directory_source}/BDI_PGC3_EUR.sumstats.gz;\
  ${path_directory_source}/EPILEPSYgen_ILAEC_TE.sumstats.gz"
)
comparisons+=(
  "${path_directory_product}/bipolar2_pgc3_against_epilepsy;\
  ${path_directory_source}/BDII_PGC3_EUR.sumstats.gz;\
  ${path_directory_source}/EPILEPSYgen_ILAEC_TE.sumstats.gz"
)


###########################################################################
# Execute procedure.

# Organize information for batch instances.
for comparison in "${comparisons[@]}"; do
  # Separate fields from instance.
  # [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]
  IFS=";" read -r -a array <<< "${comparison}"
  path_file_base_product="${array[0]}"
  path_file_source_primary="${array[1]}"
  path_file_source_secondary="${array[2]}"
  # Estimate Genetic Correlation by LDSC.
  /usr/bin/bash "${path_file_script}" \
  $path_file_source_primary \
  $path_file_source_secondary \
  $path_file_base_product \
  $path_directory_disequilibrium \
  $threads \
  $report
  # Report.
  if [[ "$report" == "true" ]]; then
    echo "----------"
    echo "Script path:"
    echo $path_script
    echo "Primary source file path:"
    echo $path_file_source_primary
    echo "Secondary source file path:"
    echo $path_file_source_secondary
    echo "Product file path:"
    echo $path_file_base_product
    echo "----------"
  fi
done

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script complete:"
  echo "4_estimate_gwas_genetic_correlation_ldsc.sh"
  echo "----------"
fi



#
