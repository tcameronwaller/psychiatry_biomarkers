#!/bin/bash

###########################################################################
###########################################################################
###########################################################################
# TODO: TCW; 15 November 2022
# TODO: verify the file paths and file names against the "table_bipolar_body_gwas_summary.ods"

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
path_gwas_summaries_waller=$(<"./gwas_summaries_waller_metabolism.txt")
path_gwas_summaries_team=$(<"./gwas_summaries_team.txt")
path_bipolar_batzler=$(<"./mayo_bipolar_batzler.txt")
path_bipolar_coombes=$(<"./mayo_bipolar_coombes.txt")

path_dock="$path_process/dock"
path_directory_product="${path_dock}/gwas_access"

# Scripts.

# Initialize directories.
#rm -r $path_directory_product
mkdir -p $path_directory_product

###########################################################################
# Execute procedure.

# Organize multi-dimensional array of information about studies.
studies=()
# [regression type] ; [full path to source file of GWAS summary statistics] ; [full path to product file of GWAS summary statistics]

##########

studies+=(
  "${path_gwas_summaries_team}/REFORMATTED_CURRENT/BMI_GIANTUKB_EUR.txt.gz;\
  ${path_directory_product}/bmi_giant_ukb.txt.gz"
)
studies+=(
  "${path_gwas_summaries_team}/REFORMATTED_CURRENT/BMI_GIANT_EUR.txt.gz;\
  ${path_directory_product}/bmi_giant.txt.gz"
)
studies+=(
  "${path_bipolar_batzler}/UKB_BMIint/output/results/PLINK/body_log.controls.18pcSiteArrayAdj.MAF0.01.glm.linear.gz;\
  ${path_directory_product}/bmi_bipolar_control_ukb.txt.gz"
)
studies+=(
  "${path_bipolar_batzler}/UKB_BMIint/output/results/PLINK/META/body_log.cases.18pcSiteArrayAdj.MAF0.01.glm.linear.Metasoft.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_ukb_metasoft.txt.gz"
)
studies+=(
  "${path_bipolar_batzler}/UKB_BMIint/output/results/R/META/body_log.cases.18pcSiteArrayAdj.MAF0.01.glm.linear.RMA.csv.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_ukb_rma.txt.gz"
)
studies+=(
  "${path_bipolar_batzler}/UKB_BMIint/output/results/R/META/body_log.cases.18pcSiteArrayAdj.MAF0.01.glm.linear.RMAFE.csv.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_ukb_rmafe.txt.gz"
)
studies+=(
  "${path_bipolar_batzler}/PGC_Analysis/metaAnalysis_imputed/BMI.pcAdj.assoc.linear/ALL_BMI.pcAdj.assoc.linear_P_MAresults.csv.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_ma.txt.gz"
)
studies+=(
  "${path_bipolar_batzler}/PGC_Analysis/metaAnalysis_imputed/BMI.pcAdj.assoc.linear/ALL_BMI.pcAdj.assoc.linear_P_MAresultsFE.csv.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_mafe.txt.gz"
)
studies+=(
  "${path_bipolar_coombes}/PGC/PGC_BMI/FUMA/ALL_BMI.pcAdj.assoc.linear_P_MAresultsFE.FUMA.txt.gz;\
  ${path_directory_product}/bmi_bipolar_case_pgc_mafe_fuma.txt.gz"
)
studies+=(
  "${path_bipolar_batzler}/UKB_BMIint/output/results/PLINK/body_log.cases.18pcSiteArrayAdj.MAF0.01.glm.linear.gz;\
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
  # Copy the source file path to the product file path.
  cp $path_file_source $path_file_product
done

################################################################################
# Report.
if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "1_access_copy_gwas.sh"
  echo "----------"
fi



#
