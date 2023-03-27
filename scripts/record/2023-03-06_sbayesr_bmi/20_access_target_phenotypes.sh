#!/bin/bash

################################################################################
# Author: T. Cameron Waller
# Date, first execution: 27 March 2023
# Date, last execution: 27 March 2023
# Review: TCW; __ March 2023
################################################################################
# Note



################################################################################
# Organize paths.

# Directories.
cd ~/paths
path_directory_genotypes=$(<"./genotypes_mayo_bipolar_disorder_1_2_merge.txt")
path_directory_genotype_pca_all=$(<"./genotype_pca_mayo_bipolar_disorder_1_2_merge_all.txt")
path_directory_genotype_pca_europe=$(<"./genotype_pca_mayo_bipolar_disorder_1_2_merge_europe.txt")
path_directory_phenotypes=$(<"./phenotypes_mayo_bipolar_disorder.txt")
path_directory_process=$(<"./process_psychiatric_metabolism.txt")
path_directory_dock="${path_directory_process}/dock" # parent directory for procedural reads and writes
path_directory_product="${path_directory_dock}/phenotypes_mayo_bipolar_disorder_1_2_merge"

# Files.
name_file_identifier="210421_id_matching_gwas.csv" # file date: 21 April 2021
name_file_genetic_sex_case="MERGED.maf0.01.dosR20.8.noDups.fam" # file date: 20 May 2022
name_file_genotype_pca="Top20_PCs.csv" # file date: 23 May 2022
name_file_genotype_pca_all="Top20_PCs_all.csv" # file date: 23 May 2022
name_file_genotype_pca_europe="Top20_PCs_europe.csv" # file date: 23 May 2022

name_file_phenotype="220513_BP_phenotypes.csv" # file date: 13 May 2022; Richard S. Pendegraft prepared and shared this file on 13 May 2022
name_file_phenotype_supplement_1="220325_BP_phenotypes.csv" # file date: 31 March 2021
name_file_phenotype_supplement_2="211221_BP_phenotypes.csv" # file date: 23 December 2021
name_file_phenotype_supplement_3="210902_BP_phenotypes.csv" # file date: 2 September 2021
name_file_phenotype_supplement_4="210609_BP_phenotypes.csv" # file date: 9 June 2021
name_file_phenotype_supplement_5="210422_BP_phenotypes.csv" # file date: 5 May 2021
name_file_phenotype_supplement_6="210330_BP_phenotypes.csv" # file date: 13 April 2021
name_file_phenotype_supplement_7="phen_bp.csv" # file date: 5 September 2022

path_file_identifier_source="${path_directory_phenotypes}/${name_file_identifier}"
path_file_genetic_sex_case_source="${path_directory_genotypes}/${name_file_genetic_sex_case}"
path_file_genotype_pca_all_source="${path_directory_genotype_pca_all}/${name_file_genotype_pca}"
path_file_genotype_pca_europe_source="${path_directory_genotype_pca_europe}/${name_file_genotype_pca}"
path_file_phenotype_source="${path_directory_phenotypes}/${name_file_phenotype}"
path_file_phenotype_supplement_1_source="${path_directory_phenotypes}/${name_file_phenotype_supplement_1}"
path_file_phenotype_supplement_2_source="${path_directory_phenotypes}/${name_file_phenotype_supplement_2}"
path_file_phenotype_supplement_3_source="${path_directory_phenotypes}/${name_file_phenotype_supplement_3}"
path_file_phenotype_supplement_4_source="${path_directory_phenotypes}/${name_file_phenotype_supplement_4}"
path_file_phenotype_supplement_5_source="${path_directory_phenotypes}/${name_file_phenotype_supplement_5}"
path_file_phenotype_supplement_6_source="${path_directory_phenotypes}/${name_file_phenotype_supplement_6}"
path_file_phenotype_supplement_7_source="${path_directory_phenotypes}/${name_file_phenotype_supplement_7}"



path_file_identifier_product="${path_directory_product}/${name_file_identifier}"
path_file_phenotype_product="${path_directory_product}/${name_file_phenotype}"
path_file_genetic_sex_case_product="${path_directory_product}/${name_file_genetic_sex_case}"
path_file_genotype_pca_all_product="${path_directory_product}/${name_file_genotype_pca_all}"
path_file_genotype_pca_europe_product="${path_directory_product}/${name_file_genotype_pca_europe}"



name_file_genotypes_prefix="MERGED.maf0.dosR20.3.noDups.chr"
name_file_genotypes_suffix=".dose.vcf.gz"
name_file_genotypes_suffix_index=".dose.vcf.gz.tbi"
#chromosome="10"
#path_file_genotypes_source="${path_directory_source}/${name_file_genotypes_prefix}${chromosome}${name_file_genotypes_suffix}"
#path_file_genotypes_product="${path_directory_product}/${name_file_genotypes_prefix}${chromosome}${name_file_genotypes_suffix}"

# Initialize directories.
rm -r $path_directory_product
mkdir -p $path_directory_product
cd $path_directory_product

###########################################################################
# Organize parameters.



###########################################################################
# Execute procedure.

# Echo each command to console.
set -x

cp ${path_directory_source}/${name_file_genotypes_prefix}*${name_file_genotypes_suffix}* $path_directory_product
#cp "$path_file_genotypes_source" "$path_file_genotypes_product"
#mv "$path_target_container/MERGED" "$path_target_mayo_bipolar_genotype"



#
