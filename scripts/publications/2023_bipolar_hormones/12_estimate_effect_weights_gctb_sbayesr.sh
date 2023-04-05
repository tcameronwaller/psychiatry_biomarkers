#!/bin/bash

################################################################################
################################################################################
################################################################################
# Author: T. Cameron Waller
# Date, first execution: __ February 2023
# Date, last execution: __ February 2023
################################################################################
################################################################################
################################################################################
# Note



################################################################################
################################################################################
################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_tools=$(<"./waller_tools.txt")
#path_plink2=$(<"./tools_plink2.txt")
path_gctb=$(<"./tools_gctb.txt")
path_directory_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_source="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/gwas_extra_process"
path_directory_product="${path_directory_dock}/hormone_genetics_tcw_2023-02-24/effect_weights_sbayesr"
path_directory_ld_matrix="${path_directory_product}/ukbEURu_hm3_shrunk_sparse"

# Files.
#path_file_gwas_source="${path_directory_source}/32042192_ruth_2020_testosterone_female.txt.gz"
#path_file_gwas_product="${path_directory_product}/32042192_ruth_2020_testosterone_female.ma"
path_file_gwas_source="${path_directory_source}/32769997_zhou_2020_tsh.txt.gz"
path_file_gwas_product="${path_directory_product}/32769997_zhou_2020_tsh.ma"
path_file_ld_matrix_source="${path_directory_reference}/gctb/ukbEURu_hm3_sparse.zip"
path_file_ld_matrix_product="${path_directory_product}/ukbEURu_hm3_sparse.zip"
name_file_ld_matrix_prefix="ukbEURu_hm3_chr"
name_file_ld_matrix_suffix="_v3_50k.ldm.sparse"
name_file_product_prefix="32769997_zhou_2020_tsh_"
name_file_product_suffix="_2023-02-27" # Suffix must not be an empty string.

# Scripts.
path_script_gwas_format="${path_directory_process}/promiscuity/scripts/gctb/constrain_translate_gwas_standard_to_gctb.sh"
path_script_submit_batch="${path_directory_process}/promiscuity/scripts/gctb/1_submit_batch_gctb_sbayesr_chromosomes.sh"
path_script_batch_run_sbayesr="${path_directory_process}/promiscuity/scripts/gctb/2_run_batch_gctb_sbayesr.sh"
path_script_run_sbayesr="${path_directory_process}/promiscuity/scripts/gctb/run_gctb_sbayesr.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product



###########################################################################
# Organize parameters.

chromosome_x="false"
threads=1
report="true"


# TODO: TCW; 27 February 2023
# TODO: need a new parameter for each GWAS sum stats <-- "observations_variant"
# TODO: derive from the "fill_observations" parameter in the GWAS format translation parameter table.



################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo "call_submit_batch_gctb_sbayesr.sh"
  echo "----------"
fi


###########################################################################
# Execute procedure.

##########
# 1. Prepare GWAS summary statistics.

if true; then
  /usr/bin/bash $path_script_gwas_format \
  $path_file_gwas_source \
  $path_file_gwas_product \
  $report
fi



##########
# 2. Prepare LD matrices.

if true; then
  cp $path_file_ld_matrix_source $path_file_ld_matrix_product
  unzip $path_file_ld_matrix_product -d $path_directory_product
fi



##########
# 3. Prepare and submit batch of jobs for processing on each chromosome.

if true; then
  /usr/bin/bash $path_script_submit_batch \
  $path_file_gwas_product \
  $path_directory_ld_matrix \
  $name_file_ld_matrix_prefix \
  $name_file_ld_matrix_suffix \
  $path_directory_product \
  $name_file_product_prefix \
  $name_file_product_suffix \
  $chromosome_x \
  $path_script_batch_run_sbayesr \
  $path_script_run_sbayesr \
  $path_gctb \
  $threads \
  $report
fi



##########
# 4. Remove temporary files.

# rm -r $path_directory_ld_matrix



#
