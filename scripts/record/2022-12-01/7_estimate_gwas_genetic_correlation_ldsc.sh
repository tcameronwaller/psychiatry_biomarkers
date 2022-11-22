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
path_process=$(<"./process_psychiatric_metabolism.txt")

path_dock="$path_process/dock"
path_directory_source="${path_dock}/bipolar_body/gwas_munge_ldsc"
path_directory_product="${path_dock}/bipolar_body/gwas_genetic_correlation_ldsc"
path_directory_reference="${path_dock}/bipolar_body/reference_ldsc"
path_directory_disequilibrium="${path_directory_reference}/disequilibrium/eur_w_ld_chr"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_directory_ldsc="${path_promiscuity_scripts}/utility/ldsc"
path_script="${path_directory_ldsc}/estimate_gwas_genetic_correlation_ldsc.sh"

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
  "${path_directory_product}/bmi_pgc_ukb_bipolar_case_metasoft_versus_bmi_giant;\
  ${path_directory_source}/bmi_bipolar_case_pgc_ukb_metasoft.sumstats.gz;\
  ${path_directory_source}/bmi_giant.sumstats.gz"
)
comparisons+=(
  "${path_directory_product}/bmi_pgc_ukb_bipolar_case_rma_versus_bmi_giant;\
  ${path_directory_source}/bmi_bipolar_case_pgc_ukb_rma.sumstats.gz;\
  ${path_directory_source}/bmi_giant.sumstats.gz"
)
comparisons+=(
  "${path_directory_product}/bmi_pgc_ukb_bipolar_case_rmafe_versus_bmi_giant;\
  ${path_directory_source}/bmi_bipolar_case_pgc_ukb_rmafe.sumstats.gz;\
  ${path_directory_source}/bmi_giant.sumstats.gz"
)
comparisons+=(
  "${path_directory_product}/bmi_pgc_ukb_bipolar_case_rmafe_versus_bmi_ukb_bipolar_control;\
  ${path_directory_source}/bmi_bipolar_case_pgc_ukb_rmafe.sumstats.gz;\
  ${path_directory_source}/bmi_bipolar_control_ukb.sumstats.gz"
)
comparisons+=(
  "${path_directory_product}/bmi_pgc_bipolar_case_mafe_fuma_versus_bmi_giant_ukb;\
  ${path_directory_source}/bmi_bipolar_case_pgc_mafe_fuma.sumstats.gz;\
  ${path_directory_source}/bmi_giant_ukb.sumstats.gz"
)
comparisons+=(
  "${path_directory_product}/bmi_pgc_bipolar_case_mafe_fuma_versus_bmi_ukb_bipolar_control;\
  ${path_directory_source}/bmi_bipolar_case_pgc_mafe_fuma.sumstats.gz;\
  ${path_directory_source}/bmi_bipolar_control_ukb.sumstats.gz"
)


###########################################################################
# Execute procedure.

# Organize information for batch instances.
for study_details in "${studies[@]}"; do
  # Separate fields from instance.
  # [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]
  IFS=";" read -r -a array <<< "${study_details}"
  path_file_base_product="${array[0]}"
  path_file_source_primary="${array[1]}"
  path_file_source_secondary="${array[2]}"
  # Estimate Genetic Correlation by LDSC.
  /usr/bin/bash "${path_script}" \
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
  echo "7_estimate_gwas_genetic_correlation_ldsc.sh"
  echo "----------"
fi



#
