#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# Review: TCW; __ November 2022

# next script in pipeline
# 1. translate format of all GWAS summary stats files to the common team format
# 2. translate format of all GWAS summary stats files to the LDSC format

###########################################################################
###########################################################################
###########################################################################



################################################################################
# Organize paths.
# Read private, local file paths.
cd ~/paths
path_process=$(<"./process_psychiatric_metabolism.txt")

path_dock="$path_process/dock"
path_directory_source="${path_dock}/bipolar_body/gwas_format_ldsc"
path_directory_product="${path_dock}/bipolar_body/gwas_munge_ldsc"
path_directory_reference="${path_dock}/bipolar_body/reference_ldsc"
path_directory_disequilibrium="${path_directory_reference}/disequilibrium/eur_w_ld_chr"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_directory_ldsc="${path_promiscuity_scripts}/utility/ldsc"
path_script="${path_directory_ldsc}/munge_gwas_ldsc.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product

# Report.
report="true"

###########################################################################
# Execute procedure.

##########
# Organize multi-dimensional array of information about studies.
studies=()
# [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]

studies+=(
  "${path_directory_source}/bmi_giant_ukb.txt.gz;\
  ${path_directory_product}/bmi_giant_ukb.txt.gz"
)
studies+=(
  "${path_directory_source}/bmi_giant.txt.gz;\
  ${path_directory_product}/bmi_giant.txt.gz"
)
studies+=(
  "${path_directory_source}/bmi_bipolar_control_ukb.txt.gz;\
  ${path_directory_product}/bmi_bipolar_control_ukb.txt.gz"
)
studies+=(
  "${path_directory_source}/bmi_bipolar_case_pgc_ukb_metasoft.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_ukb_metasoft.txt.gz"
)
studies+=(
  "${path_directory_source}/bmi_bipolar_case_pgc_ukb_rma.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_ukb_rma.txt.gz"
)
studies+=(
  "${path_directory_source}/bmi_bipolar_case_pgc_ukb_rmafe.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_ukb_rmafe.txt.gz"
)
studies+=(
  "${path_directory_source}/bmi_bipolar_case_pgc_ma.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_ma.txt.gz"
)
studies+=(
  "${path_directory_source}/bmi_bipolar_case_pgc_mafe.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_mafe.txt.gz"
)
studies+=(
  "${path_directory_source}/bmi_bipolar_case_pgc_mafe_fuma.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_mafe_fuma.txt.gz"
)
studies+=(
  "${path_directory_source}/bmi_bipolar_case_ukb.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_ukb.txt.gz"
)

################################################################################

# Organize information for batch instances.
for study_details in "${studies[@]}"; do
  # Separate fields from instance.
  # [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]
  IFS=";" read -r -a array <<< "${study_details}"
  path_file_source="${array[0]}"
  path_file_product="${array[1]}"
  # Translate GWAS summary statistics to team standard format.
  /usr/bin/bash "${path_script}" \
  $path_file_source \
  $path_file_product \
  $report
  # Report.
  if [[ "$report" == "true" ]]; then
    echo "----------"
    echo "Script path:"
    echo $path_script
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
  echo "3_format_gwas_ldsc.sh"
  echo "----------"
fi



#
