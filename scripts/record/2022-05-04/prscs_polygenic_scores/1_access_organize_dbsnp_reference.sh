#!/bin/bash

###########################################################################
# Organize paths.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_bcftools=$(<"./tools_bcftools.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="${path_process}/dock"
path_parameters="${path_dock}/parameters"
path_translations_chromosomes_dbsnp="${path_parameters}/promiscuity/translations_chromosomes_refseq_grch38p14.txt"
path_dbsnp_raw_container="${path_dock}/access/dbsnp_reference_raw"
#path_dbsnp_raw_file="${path_dock}/access/dbsnp_reference_raw/GCF_000001405.39.gz"
path_dbsnp_format_container="${path_dock}/access/dbsnp_reference_format_thread_12"
#path_dbsnp_format_file="${path_dock}/access/dbsnp_reference_format/GCF_000001405.39.gz"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_access_dbsnp="${path_promiscuity_scripts}/utility/access_dbsnp_reference.sh"
path_script_chromosome_in_vcf="${path_promiscuity_scripts}/utility/bcftools/translate_chromosomes_in_vcf.sh"

###########################################################################
# Execute procedure.
###########################################################################

###########################################################################
# Access raw files from dbSNP.

if false; then
  /usr/bin/bash "${path_script_access_dbsnp}" \
  $path_dbsnp_raw_container
fi

###########################################################################
# Format genotype information in VCF.

if true; then
  # Initialize directory.
  rm -r $path_dbsnp_format_container
  mkdir -p $path_dbsnp_format_container
  # Organize specific paths and parameters.
  path_vcf_source="${path_dbsnp_raw_container}/GCF_000001405.39.gz"
  path_vcf_product="${path_dbsnp_format_container}/GCF_000001405.39.gz" # determine suffix from BCFTools argument
  report="true"
  # Convert information from genotype files in VCF format to BIM format.
  /usr/bin/bash "${path_script_chromosome_in_vcf}" \
  $path_translations_chromosomes_dbsnp \
  $path_vcf_source \
  $path_vcf_product \
  $path_bcftools \
  $report
fi
