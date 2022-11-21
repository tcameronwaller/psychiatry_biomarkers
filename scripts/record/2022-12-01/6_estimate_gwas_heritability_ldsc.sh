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
path_directory_product="${path_dock}/bipolar_body/gwas_heritability_ldsc"
path_directory_reference="${path_dock}/bipolar_body/reference_ldsc"
path_directory_disequilibrium="${path_directory_reference}/disequilibrium/eur_w_ld_chr"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_directory_ldsc="${path_promiscuity_scripts}/utility/ldsc"
path_script="${path_directory_ldsc}/estimate_gwas_heritability_ldsc.sh"

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
studies=()
# [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]

studies+=(
  "${path_directory_source}/bmi_giant_ukb.sumstats.gz;\
  ${path_directory_product}/bmi_giant_ukb"
)
studies+=(
  "${path_directory_source}/bmi_giant.sumstats.gz;\
  ${path_directory_product}/bmi_giant"
)
studies+=(
  "${path_directory_source}/bmi_bipolar_control_ukb.sumstats.gz;\
  ${path_directory_product}/bmi_bipolar_control_ukb"
)
studies+=(
  "${path_directory_source}/bmi_bipolar_case_pgc_ukb_metasoft.sumstats.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_ukb_metasoft"
)
studies+=(
  "${path_directory_source}/bmi_bipolar_case_pgc_ukb_rma.sumstats.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_ukb_rma"
)
studies+=(
  "${path_directory_source}/bmi_bipolar_case_pgc_ukb_rmafe.sumstats.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_ukb_rmafe"
)
if false; then
  # Even after extensive trouble-shooting, LDSC Munge is unable to process these
  # GWAS summary statistics.
  # TCW; 21 November 2022
  studies+=(
    "${path_directory_source}/bmi_bipolar_case_pgc_ma.sumstats.gz;\
    ${path_directory_product}/bmi_bipolar_case_pgc_ma"
  )
  studies+=(
    "${path_directory_source}/bmi_bipolar_case_pgc_mafe.sumstats.gz;\
    ${path_directory_product}/bmi_bipolar_case_pgc_mafe"
  )
fi
studies+=(
  "${path_directory_source}/bmi_bipolar_case_pgc_mafe_fuma.sumstats.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_mafe_fuma"
)
studies+=(
  "${path_directory_source}/bmi_bipolar_case_ukb.sumstats.gz;\
  ${path_directory_product}/bmi_bipolar_case_ukb"
)

###########################################################################
# Execute procedure.

# Organize information for batch instances.
for study_details in "${studies[@]}"; do
  # Separate fields from instance.
  # [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]
  IFS=";" read -r -a array <<< "${study_details}"
  path_file_source="${array[0]}"
  path_file_base_product="${array[1]}"
  # Translate GWAS summary statistics to team standard format.
  /usr/bin/bash "${path_script}" \
  $path_file_source \
  $path_file_base_product \
  $path_directory_disequilibrium \
  $threads \
  $report
  # Report.
  if [[ "$report" == "true" ]]; then
    echo "----------"
    echo "Script path:"
    echo $path_script
    echo "Source file path:"
    echo $path_file_source
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
  echo "6_estimate_gwas_heritability_ldsc.sh"
  echo "----------"
fi



#
