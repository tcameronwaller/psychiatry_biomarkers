#!/bin/bash

################################################################################
################################################################################
################################################################################
# Author: T. Cameron Waller
# Date, first execution: 1 March 2023
# Date, last execution: 1 March 2023
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
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock"
#path_directory_parameters="${path_directory_dock}/parameters/psychiatric_metabolism"
path_directory_source="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-01/gwas_format_standard"
path_directory_product="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-01/gwas_vcf_process"

# Files.
path_file_batch_instances="${path_directory_product}/batch_instances.txt"

# Scripts.
#path_promiscuity_scripts="${path_directory_process}/promiscuity/scripts"
#path_script_submit_batch="${path_promiscuity_scripts}/gwas_clean/1_submit_batch_pipe_gwas_clean.sh"
path_file_script_run_batch="${path_directory_process}/promiscuity/scripts/gwas_clean/2_run_batch_pipe_gwas_clean.sh"
path_file_script_pipe_gwas_clean="${path_directory_process}/promiscuity/scripts/gwas_clean/pipe_gwas_clean.sh"

# Initialize files.
rm $path_file_batch_instances

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product

################################################################################
# Organize parameters.

threads=1
report="true"



################################################################################
# Execute procedure.

# Define variables.
name="BMI_GIANTUKB_EUR"
type="linear"
# Organize paths.
path_file_gwas_source="${path_directory_source}/${name}.txt.gz"
path_file_gwas_product="${path_directory_product}/${name}.txt.gz"
# Define and append a new batch instance.
instance="${path_file_gwas_source};${path_file_gwas_product};${type}"
echo $instance >> $path_file_batch_instances



################################################################################
# Submit batch instances to cluster scheduler.

# Read batch instances.
readarray -t batch_instances < $path_file_batch_instances
batch_instances_count=${#batch_instances[@]}
echo "----------"
echo "count of batch instances: " $batch_instances_count
echo "first batch instance: " ${batch_instances[0]} # notice base-zero indexing
echo "last batch instance: " ${batch_instances[$batch_instances_count - 1]}

# Execute batch with grid scheduler.
if true; then
  # Submit array batch to Sun Grid Engine.
  # Array batch indices must start at one (not zero).
  qsub -t 1-${batch_instances_count}:1 \
  -o "${path_directory_product}/batch_out.txt" \
  -e "${path_directory_product}/batch_error.txt" \
  $path_file_script_run_batch \
  $path_file_batch_instances \
  $batch_instances_count \
  $path_file_script_pipe_gwas_clean \
  $threads \
  $report
fi



#
