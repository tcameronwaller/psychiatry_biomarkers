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
path_directory_source="${path_dock}/bipolar_body/gwas_access"
path_directory_product="${path_dock}/bipolar_body/gwas_format_team"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_directory_format_team="${path_promiscuity_scripts}/gwas_process/format_gwas_team"

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

# format group: already in team standard format
# delimiter: white space
# SNP CHR BP A1 A2 A1AF BETA SE P N Z INFO NCASE NCONT
studies+=(
  "none;\
  ${path_directory_source}/bmi_giant_ukb.txt.gz;\
  ${path_directory_product}/bmi_giant_ukb.txt.gz"
)
studies+=(
  "none;\
  ${path_directory_source}/bmi_giant.txt.gz;\
  ${path_directory_product}/bmi_giant.txt.gz"
)

# format group: PLINK2 linear
# delimiter: white space
# CHROM POS ID REF ALT A1 A1_FREQ TEST OBS_CT BETA SE T_STAT P
studies+=(
  "${path_directory_format_team}/translate_gwas_plink_linear.sh;\
  ${path_directory_source}/bmi_bipolar_control_ukb.txt.gz;\
  ${path_directory_product}/bmi_bipolar_control_ukb.txt.gz"
)
studies+=(
  "${path_directory_format_team}/translate_gwas_plink_linear.sh;\
  ${path_directory_source}/bmi_bipolar_case_ukb.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_ukb.txt.gz"
)

# format group: Metasoft
# delimiter: white space; tab
# RSID CHROM POS A1 A2 inMeta #STUDY PVALUE_FE BETA_FE STD_FE PVALUE_RE BETA_RE
# STD_RE PVALUE_RE2 STAT1_RE2 STAT2_RE2 PVALUE_BE I_SQUARE Q PVALUE_Q TAU_SQUARE
# PVALUES_OF_STUDIES(Tab_delimitered) MVALUES_OF_STUDIES(Tab_delimitered)
studies+=(
  "${path_directory_format_team}/translate_gwas_metasoft_body_bipolar.sh;\
  ${path_directory_source}/bmi_bipolar_case_pgc_ukb_metasoft.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_ukb_metasoft.txt.gz"
)

# format group: R rma
# delimiter: comma
# RSID CHROM POS A1 A2 inMeta QEp b se pval
studies+=(
  "${path_directory_format_team}/translate_gwas_r_metafor_body_bipolar.sh;\
  ${path_directory_source}/bmi_bipolar_case_pgc_ukb_rma.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_ukb_rma.txt.gz"
)
studies+=(
  "${path_directory_format_team}/translate_gwas_r_metafor_body_bipolar.sh;\
  ${path_directory_source}/bmi_bipolar_case_pgc_ukb_rmafe.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_ukb_rmafe.txt.gz"
)

# format group: PGC
# delimiter: comma
# SNP A1 A2 QEp b se pval
studies+=(
  "${path_directory_format_team}/translate_gwas_pgc_body_bipolar.sh;\
  ${path_directory_source}/bmi_bipolar_case_pgc_ma.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_ma.txt.gz"
)
studies+=(
  "${path_directory_format_team}/translate_gwas_pgc_body_bipolar.sh;\
  ${path_directory_source}/bmi_bipolar_case_pgc_mafe.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_mafe.txt.gz"
)

# format group: FUMA
# delimiter: white space
# SNP CHR POS A1 A2 QEp b se pval
studies+=(
  "${path_directory_format_team}/translate_gwas_pgc_fuma_body_bipolar.sh;\
  ${path_directory_source}/bmi_bipolar_case_pgc_mafe_fuma.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_mafe_fuma.txt.gz"
)

################################################################################

# Organize information for batch instances.
for study_details in "${studies[@]}"; do
  # Separate fields from instance.
  # [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]
  IFS=";" read -r -a array <<< "${study_details}"
  path_script="${array[0]}"
  path_file_source="${array[1]}"
  path_file_product="${array[2]}"
  if [[ "$path_script" != "none" ]]; then
    # Translate GWAS summary statistics to team standard format.
    /usr/bin/bash "${path_script}" \
    $path_file_source \
    $path_file_product \
    $report
  else
    # Copy the GWAS summary statistics alread in team standard format.
    cp $path_file_source $path_file_product
  fi
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
  echo "2_format_standard_gwas.sh"
  echo "----------"
fi



#
