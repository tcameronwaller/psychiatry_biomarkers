#!/bin/bash

################################################################################
# Organize argument variables.

path_table_phenotypes_covariates=${1} # full path to file for table with phenotypes and covariates
path_report=${2} # full path to parent directory for GWAS summary statistics
analysis=${3} # unique name for association analysis
phenotypes=${4} # names of table's column or columns for single or multiple phenotypes, dependent variables
covariates=${5} # name of table's columns for covariates, independent variables
threads=${6} # count of processing threads to use
maf=${7} # minor allele frequency threshold filter
chromosomes=${8} # count of sequential chromosomes

###########################################################################
# Organize variables.

# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_plink2=$(<"./tools_plink2.txt")
path_ukb_genotype=$(<"./ukbiobank_genotype.txt")

# Iterate on chromosomes.
start=1
end=$chromosomes
for (( index=$start; index<=$end; index++ ))
do

  #echo "chromosome: ${index}"
  # Set directory.
  path_chromosome="$path_report/chromosome_${index}"
  # Determine whether the temporary directory structure already exists.
  if [ ! -d $path_chromosome ]; then
      # Directory does not already exist.
      # Create directory.
      mkdir -p $path_chromosome
  fi
  cd $path_chromosome

  if true; then
    # Call PLINK2.
    # 90,000 Mebibytes (MiB) is 94.372 Gigabytes (GB)
    # --pfilter drops SNPs with null p-values and any beyond threshold (such as 1)
    $path_plink2 \
    --memory 90000 \
    --threads $threads \
    --bgen $path_ukb_genotype/Chromosome/ukb_imp_chr${index}_v3.bgen \
    --sample $path_ukb_genotype/Chromosome/ukb46237_imp_chr${index}_v3_s487320.sample \
    --keep $path_table_phenotypes_covariates \
    --maf $maf \
    --freq --glm hide-covar \
    --pfilter 1 \
    --pheno $path_table_phenotypes_covariates \
    --pheno-name $phenotypes \
    --covar $path_table_phenotypes_covariates \
    --covar-name $covariates \
    --out report
  fi

done
