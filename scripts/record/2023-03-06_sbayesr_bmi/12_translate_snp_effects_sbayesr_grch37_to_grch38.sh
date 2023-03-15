#!/bin/bash

################################################################################
################################################################################
################################################################################
# Author: T. Cameron Waller
# Date, first execution: 10 March 2023
# Date, last execution: 15 March 2023
################################################################################
################################################################################
################################################################################
# Note


# TODO: TCW; 10 March 2023
# 0. start using a more legit reference genome sequence in GWAS2VCF
# 1. translate to BED format for CrossMap
# 2. map in CrossMap GRCh37 to GRCh38
# 3. translate to format for PLINK2 --score <-- Might not need to do this... just tell PLINK which columns?
# 4. apply PLINK2 --score



################################################################################
################################################################################
################################################################################



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes

path_directory_source="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-01/combination_effect_weights_sbayesr"
path_directory_product="${path_directory_dock}/test_sbayesr_body_mass_tcw_2023-03-01/sbayesr_snp_effects_grch38"

# Files.
path_file_grch37="${path_directory_source}/BMI_GIANTUKB_EUR_tcw_2023-03-01_chromosomes.snpRes"
path_file_grch37_ucsc_bed="${path_directory_product}/BMI_GIANTUKB_EUR_grch37.bed.gz"
#path_file_chain_grch37_to_grch38="${path_directory_reference}/assembly_chains/ucsc/hg19ToHg38.over.chain.gz"
path_file_chain_grch37_to_grch38="${path_directory_reference}/crossmap/ensembl/GRCh37_to_GRCh38.chain.gz"
path_file_grch38_ucsc_bed="${path_directory_product}/BMI_GIANTUKB_EUR_grch38.bed.gz"
path_file_grch38="${path_directory_product}/BMI_GIANTUKB_EUR_grch38_standard.txt.gz"

# Scripts.
path_script_ucsc_bed="${path_directory_process}/promiscuity/scripts/gctb/translate_snp_effects_sbayesr_to_ucsc_bed.sh"
path_script_map="${path_directory_process}/promiscuity/scripts/crossmap/map_genomic_feature_bed.sh"
path_script_standard="${path_directory_process}/promiscuity/scripts/gctb/translate_snp_effects_ucsc_bed_to_standard.sh"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product



###########################################################################
# Organize parameters.

threads=1
report="true"

###########################################################################
# Execute procedure.

##########
# Translate SBayesR SNP effect weights to format for CrossMap.

if true; then
  /usr/bin/bash $path_script_ucsc_bed \
  $path_file_grch37 \
  $path_file_grch37_ucsc_bed \
  $report
fi

##########
# Translate SBayesR SNP effect weights in CrossMap from GRCh37 to GRCh38.

if true; then
  /usr/bin/bash $path_script_map \
  $path_file_grch37_ucsc_bed \
  $path_file_grch38_ucsc_bed \
  $path_file_chain_grch37_to_grch38 \
  $threads \
  $report
fi

##########
# Translate SBayesR SNP effect weights in CrossMap from GRCh37 to GRCh38.

if true; then
  /usr/bin/bash $path_script_standard \
  $path_file_grch38_ucsc_bed \
  $path_file_grch38 \
  $report
fi


################################################################################
# Report.

if [[ "$report" == "true" ]]; then
  echo "----------"
  echo "Script:"
  echo "12_translate_snp_effects_sbayesr_grch37_to_grch38.sh"
  echo "----------"
fi



#
