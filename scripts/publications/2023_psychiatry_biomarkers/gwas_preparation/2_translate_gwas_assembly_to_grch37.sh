#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 24 May 2023
# Date, last execution: 19 December 2023
# Date, review: 19 December 2023
################################################################################
# Note

# The purpose of this script is to translate to human genome assembly GRCh37
# (team standard) any GWAS summary statistics that use coordinates from any
# other human genome assembly (NCBI36-hg18, GRCh38, etc).

# On 24 May 2023, TCW chose to prioritize the chain map coordinates from
# Ensembl.

# After running this script, check the product directory to make sure that this
# procedure wrote the appropriate files again, later than all others.

# Note: TCW; 13 December 2023
# 1. Description of situation
# From the publication article and supplement, it is not definitively clear
# whether the original GWAS summary statistics from study
# "37872160_williams_2023" used coordinates from human genome assembly GRCh37 or
# GRCh38. This set of GWAS summary statistics had 57,524,163 SNPs after format
# translation to standard.
# 2. Without any attempt to translate from GRCh38 to GRCh37
# On 13 December 2023 and again on 14 December 2023, without any prior attempt
# to translate from GRCh38 to GRCh37 but with prior filters on SNPs (script
# file: "filter_constrain_gwas_summary_values.sh"; 57,523,911 SNPs), the
# filtered GWAS summary statistics completed a procedure
# (script file: "fill_reference_snp_cluster_identifier.sh") that filled SNP
# rsIDs from dbSNP (GRCh37). Of the 57,523,911 SNPs originally in the GWAS
# summary statistics, 36,107,225 SNPs (62.8%) matched the dbSNP reference
# successfully.
# 3. With attempt to translate from GRCh38 to GRCh37
# On 30 November 2023, the source GWAS summary statistics with 57,524,163
# original SNPs completed an attempt to map from GRCh38 to GRCh37 (Crossmap;
# Ensembl chain) with 54,848,929 final SNPs, and these same GWAS summary
# statistics subsequently completed the GWAS2VCF procedure with 10,152,437 SNPs
# on 4 December 2023.
# Separately, on 13 December 2023, the source GWAS summary statistics with
# 57,524,163 original SNPs completed an attempt to map from GRCh38 to GRCh37
# (Crossmap; Ensembl chain) with 54,848,929 final SNPs. These same GWAS summary
# statistics (54,848,929 SNPs) then completed a procedure
# (script file: "fill_reference_snp_cluster_identifier.sh") that filled SNP
# rsIDs from dbSNP (GRCh37). Of the 54,848,929 SNPs originally in the GWAS
# summary statistics, 4,073,715 SNPs (7.4%) matched the dbSNP reference
# successfully. The fact that these SNPs matched on the basis of chromosome,
# position, reference allele, and alternate allele does not necessarily mean
# that the human genome assemblies were both correct and the same.
# 4. Assessment
# Without an attempt to map from GRCh38 to GRCh37, a much greater proportion of
# the SNPs from the original GWAS summary statistics (63% compared to 7%)
# matched successfully with SNPs from the dbSNP reference (GRCh37). This
# evidence suggests that the original GWAS summary statistics from study
# "37872160_williams_2023" used the GRCh37 assembly of the human genome.



################################################################################
# Organize paths.

# Identifiers or designators of parameter version and preparation batch.
identifier_preparation="gwas_2023-12-19_alcohol_sex_test_ldsc_2023-12-19"
identifier_parameter="tcw_2023-12-19_alcohol_sex_test"

# Directories.
cd ~/paths
path_directory_reference=$(<"./reference_tcw.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_source="${path_directory_dock}/${identifier_preparation}/1_gwas_format_standard"
path_directory_product="${path_directory_dock}/${identifier_preparation}/2_gwas_assembly_grch37"

# Files.
#path_file_chain_ncbi36_to_grch37="${path_directory_reference}/crossmap/ucsc/hg18ToHg19.over.chain.gz"
#path_file_chain_grch37_to_grch38="${path_directory_reference}/crossmap/ucsc/hg19ToHg38.over.chain.gz"
#path_file_chain_grch38_to_grch37="${path_directory_reference}/crossmap/ucsc/hg38ToHg19.over.chain.gz"
path_file_chain_ncbi36_to_grch37="${path_directory_reference}/crossmap/ensembl/NCBI36_to_GRCh37.chain.gz"
path_file_chain_grch37_to_grch38="${path_directory_reference}/crossmap/ensembl/GRCh37_to_GRCh38.chain.gz"
path_file_chain_grch38_to_grch37="${path_directory_reference}/crossmap/ensembl/GRCh38_to_GRCh37.chain.gz"

# Scripts.
path_directory_partner_scripts="${path_directory_process}/partner/scripts"
path_script_map_assembly="${path_directory_partner_scripts}/crossmap/map_gwas_standard_format_bed.sh"

# Initialize directories.
rm -r $path_directory_product # caution
mkdir -p $path_directory_product
cd $path_directory_product

###########################################################################
# Organize parameters.

threads=1
report="true"

################################################################################
# Execute procedure.

##########
# Copy the GWAS summary statistics from the previous process.
# Most sets of GWAS summary statistics do not need extra processing.
# Subsequent processes on a few studies will replace the appropriate files.
cp $path_directory_source/*.txt.gz $path_directory_product

##########
# Translate genomic assemblies to GRCh37.



##########
# 36635386_chen_2023 (GRCh38 to GRCh37)

if false; then

  # UCSC: 15,428,167 lines to 15,396,790 lines (TCW; 24 May 2023)
  # Ensembl: 15,428,167 lines to 15,384,987 lines (TCW; 26 September 2023)
  /usr/bin/bash $path_script_map_assembly \
  "${path_directory_source}/36635386_chen_2023_cortisol.txt.gz" \
  "${path_directory_product}/36635386_chen_2023_cortisol.txt.gz" \
  $path_file_chain_grch38_to_grch37 \
  $threads \
  $report

  # Ensembl: 15,431,007 lines to 15,387,813 lines (TCW; 26 September 2023)
  /usr/bin/bash $path_script_map_assembly \
  "${path_directory_source}/36635386_chen_2023_thyroxine_total.txt.gz" \
  "${path_directory_product}/36635386_chen_2023_thyroxine_total.txt.gz" \
  $path_file_chain_grch38_to_grch37 \
  $threads \
  $report

fi

##########
# 36477530_saunders_2022 (GRCh38 to GRCh37)

# UCSC:
# Ensembl: 13,268,541 lines to 13,254,829 lines (TCW; 31 August 2023)
/usr/bin/bash $path_script_map_assembly \
"${path_directory_source}/36477530_saunders_2022_alcohol_all.txt.gz" \
"${path_directory_product}/36477530_saunders_2022_alcohol_all.txt.gz" \
$path_file_chain_grch38_to_grch37 \
$threads \
$report

# UCSC:
# Ensembl: 13,422,054 lines to 13,408,726 lines (TCW; 31 August 2023)
/usr/bin/bash $path_script_map_assembly \
"${path_directory_source}/36477530_saunders_2022_alcohol_no_ukb.txt.gz" \
"${path_directory_product}/36477530_saunders_2022_alcohol_no_ukb.txt.gz" \
$path_file_chain_grch38_to_grch37 \
$threads \
$report

if false; then

  # UCSC:
  # Ensembl: 13,763,313 lines to 13,712,881 lines (TCW; 31 August 2023)
  /usr/bin/bash $path_script_map_assembly \
  "${path_directory_source}/36477530_saunders_2022_tobacco_all.txt.gz" \
  "${path_directory_product}/36477530_saunders_2022_tobacco_all.txt.gz" \
  $path_file_chain_grch38_to_grch37 \
  $threads \
  $report

  # UCSC:
  # Ensembl: 13,924,396 lines to 13,874,259 lines (TCW; 31 August 2023)
  /usr/bin/bash $path_script_map_assembly \
  "${path_directory_source}/36477530_saunders_2022_tobacco_no_ukb.txt.gz" \
  "${path_directory_product}/36477530_saunders_2022_tobacco_no_ukb.txt.gz" \
  $path_file_chain_grch38_to_grch37 \
  $threads \
  $report

  # UCSC:
  # Ensembl: _____ lines to ____ lines (TCW; __ November 2023)
  /usr/bin/bash $path_script_map_assembly \
  "${path_directory_source}/36477530_saunders_2022_tobacco_ever_all.txt.gz" \
  "${path_directory_product}/36477530_saunders_2022_tobacco_ever_all.txt.gz" \
  $path_file_chain_grch38_to_grch37 \
  $threads \
  $report

  # UCSC:
  # Ensembl: _____ lines to ____ lines (TCW; __ November 2023)
  /usr/bin/bash $path_script_map_assembly \
  "${path_directory_source}/36477530_saunders_2022_tobacco_ever_no_ukb.txt.gz" \
  "${path_directory_product}/36477530_saunders_2022_tobacco_ever_no_ukb.txt.gz" \
  $path_file_chain_grch38_to_grch37 \
  $threads \
  $report

  # UCSC:
  # Ensembl: _____ lines to ____ lines (TCW; __ November 2023)
  /usr/bin/bash $path_script_map_assembly \
  "${path_directory_source}/36477530_saunders_2022_tobacco_age_all.txt.gz" \
  "${path_directory_product}/36477530_saunders_2022_tobacco_age_all.txt.gz" \
  $path_file_chain_grch38_to_grch37 \
  $threads \
  $report

  # UCSC:
  # Ensembl: _____ lines to ____ lines (TCW; __ November 2023)
  /usr/bin/bash $path_script_map_assembly \
  "${path_directory_source}/36477530_saunders_2022_tobacco_age_no_ukb.txt.gz" \
  "${path_directory_product}/36477530_saunders_2022_tobacco_age_no_ukb.txt.gz" \
  $path_file_chain_grch38_to_grch37 \
  $threads \
  $report

  # UCSC:
  # Ensembl: _____ lines to ____ lines (TCW; __ November 2023)
  /usr/bin/bash $path_script_map_assembly \
  "${path_directory_source}/36477530_saunders_2022_tobacco_cessation_all.txt.gz" \
  "${path_directory_product}/36477530_saunders_2022_tobacco_cessation_all.txt.gz" \
  $path_file_chain_grch38_to_grch37 \
  $threads \
  $report

  # UCSC:
  # Ensembl: _____ lines to ____ lines (TCW; __ November 2023)
  /usr/bin/bash $path_script_map_assembly \
  "${path_directory_source}/36477530_saunders_2022_tobacco_cessation_no_ukb.txt.gz" \
  "${path_directory_product}/36477530_saunders_2022_tobacco_cessation_no_ukb.txt.gz" \
  $path_file_chain_grch38_to_grch37 \
  $threads \
  $report

fi

##########
# 34662886_backman_2021 (GRCh38 to GRCh37)

# Ensembl: 502,524 lines to 499,614 lines (TCW; 26 September 2023)
/usr/bin/bash $path_script_map_assembly \
"${path_directory_source}/34662886_backman_2021_albumin.txt.gz" \
"${path_directory_product}/34662886_backman_2021_albumin.txt.gz" \
$path_file_chain_grch38_to_grch37 \
$threads \
$report



##########
# 34017140_mbatchou_2021 (GRCh38 to GRCh37)

# Ensembl: 11,367,923 lines to 11,346,500 (TCW; 26 September 2023)
/usr/bin/bash $path_script_map_assembly \
"${path_directory_source}/34017140_mbatchou_2021_albumin.txt.gz" \
"${path_directory_product}/34017140_mbatchou_2021_albumin.txt.gz" \
$path_file_chain_grch38_to_grch37 \
$threads \
$report

if false; then

  ##########
  # 32747698_matoba_2020

  # Ensembl: _____ lines to ____ lines (TCW; __ November 2023)
  /usr/bin/bash $path_script_map_assembly \
  "${path_directory_source}/32747698_matoba_2020_europe.txt.gz" \
  "${path_directory_product}/32747698_matoba_2020_europe.txt.gz" \
  $path_file_chain_grch38_to_grch37 \
  $threads \
  $report


  ##########
  # 32581359_saevarsdottir_2020 (GRCh38 to GRCh37)

  # Ensembl: 44,690,177 lines to 43,529,207 (TCW; 26 September 2023)
  /usr/bin/bash $path_script_map_assembly \
  "${path_directory_source}/32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz" \
  "${path_directory_product}/32581359_saevarsdottir_2020_thyroid_autoimmunity.txt.gz" \
  $path_file_chain_grch38_to_grch37 \
  $threads \
  $report



  ##########
  # 24586183_medici_2014 (NCBI36-hg18 to GRCh37)

  # UCSC: 2,425,175 lines to 2,424,714 lines (TCW; 24 May 2023)
  # Ensembl: 2,425,175 lines to 2,424,987 lines (TCW; 26 September 2023)
  /usr/bin/bash $path_script_map_assembly \
  "${path_directory_source}/24586183_medici_2014_thyroid_peroxidase_antibody.txt.gz" \
  "${path_directory_product}/24586183_medici_2014_thyroid_peroxidase_antibody.txt.gz" \
  $path_file_chain_ncbi36_to_grch37 \
  $threads \
  $report

  # Ensembl: 2,425,175 lines to 2,424,987 lines (TCW; 26 September 2023)
  /usr/bin/bash $path_script_map_assembly \
  "${path_directory_source}/24586183_medici_2014_thyroid_peroxidase_reactivity.txt.gz" \
  "${path_directory_product}/24586183_medici_2014_thyroid_peroxidase_reactivity.txt.gz" \
  $path_file_chain_ncbi36_to_grch37 \
  $threads \
  $report

fi

#
