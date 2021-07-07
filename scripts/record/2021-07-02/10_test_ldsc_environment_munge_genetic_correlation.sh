
################################################################################
# This script exists to test method for genetic correlation analysis by comparison to analysis by Brandon J. Coombes, Ph.D.
################################################################################

#!/bin/bash

################################################################################
# Activate Virtual Environment.

# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_tools=$(<"./waller_tools.txt")
path_environment_ldsc="${path_tools}/python/environments/ldsc"
source "${path_environment_ldsc}/bin/activate"
which python2
sleep 5s

################################################################################
# Organize paths.

# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_ldsc=$(<"./tools_ldsc.txt")
path_gwas_summaries=$(<"./gwas_summaries_waller_metabolism.txt")
#path_gwas_summaries_team=$(<"./gwas_summaries_team.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")

# Path to accession of BMI GWAS (Yengo et al, 2018; PubMed:30124842).
path_gwas_bmi="${path_gwas_summaries}/30124842_yengo_2018/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz"

path_gwas_coombes=$(<"./gwas_pgc_bipolar_bmi_coombes.txt")
path_gwas_pgc_bipolar_bmi="${path_gwas_coombes}/ALL_BMI.pcAdj.assoc.linear_P_MAresultsFE.FUMA.txt.gz"

path_dock="${path_process}/dock"
path_gwas_ukb_bipolar_control_bmi="${path_dock}/gwas/cohorts_models/white_bipolar_disorder_control_unadjust_body_mass_index/gwas_concatenation.txt.gz"
path_gwas_ukb_bipolar_case_bmi="${path_dock}/gwas/cohorts_models/white_bipolar_disorder_case_unadjust_body_mass_index/gwas_concatenation.txt.gz"

path_test_ldsc="${path_dock}/test_ldsc_environment"
path_gwas="${path_dock}/test_ldsc_environment/gwas"
path_genetic_reference="${path_dock}/test_ldsc_environment/genetic_reference"
path_alleles="${path_genetic_reference}/alleles"
path_disequilibrium="${path_genetic_reference}/disequilibrium"

path_heritability="${path_dock}/test_ldsc_environment/heritability"
path_genetic_correlation="${path_dock}/test_ldsc_environment/genetic_correlation"

################################################################################
# Organize directories.

rm -r $path_test_ldsc
# Determine whether the temporary directory structure already exists.
if [ ! -d $path_test_ldsc ]; then
    # Directory does not already exist.
    # Create directory.
    mkdir -p $path_test_ldsc
    mkdir -p $path_gwas
    mkdir -p $path_genetic_reference
    mkdir -p $path_alleles
    mkdir -p $path_disequilibrium
    mkdir -p $heritability
    mkdir -p $path_genetic_correlation
fi

################################################################################
# Access references for LDSC.

cd $path_genetic_reference

# Definitions of Simple Nucleotide Variant alleles.
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
bunzip2 "$path_genetic_reference/w_hm3.snplist.bz2"
mv "$path_genetic_reference/w_hm3.snplist" "$path_alleles/w_hm3.snplist"
# w_hm3.snplist

# Linkage disequilibrium scores for European population.
# For simple heritability estimation.
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
tar -xjvf eur_w_ld_chr.tar.bz2 -C $path_disequilibrium
# dock/access/disequilibrium/eur_w_ld_chr/*

################################################################################
# Format GWAS summary statistics for analysis in LDSC.

cd $path_gwas

# Accession of BMI GWAS summary statistics.
echo "SNP A1 A2 N BETA P" > gwas_bmi_format.txt
zcat $path_gwas_bmi | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {print $3, toupper($4), toupper($5), $10, $7, $9}' >> gwas_bmi_format.txt
head gwas_bmi_format.txt

# BMI GWAS in Bipolar Disorder Controls in U.K. Biobank.
echo "SNP A1 A2 N BETA P" > gwas_ukb_bipolar_control_bmi_format.txt
zcat $path_gwas_ukb_bipolar_control_bmi | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {print $3, toupper($6), toupper($4), $8, $9, $12}' >> gwas_ukb_bipolar_control_bmi_format.txt
head gwas_ukb_bipolar_control_bmi_format.txt

# BMI GWAS in Bipolar Disorder Cases in U.K. Biobank.
echo "SNP A1 A2 N BETA P" > gwas_ukb_bipolar_case_bmi_format.txt
zcat $path_gwas_ukb_bipolar_case_bmi | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {print $3, toupper($6), toupper($4), $8, $9, $12}' >> gwas_ukb_bipolar_case_bmi_format.txt
head gwas_ukb_bipolar_case_bmi_format.txt

# A few SNPs do not have rsIDs and instead use chromosome and position.
# Coombes BMI in Bipolar Disorder (PGC cohort) GWAS summary statistics.
echo "SNP A1 A2 N BETA P" > gwas_pgc_bipolar_bmi_format.txt
zcat $path_gwas_pgc_bipolar_bmi | awk 'BEGIN { FS=" "; OFS=" " } NR > 1 {print $1, toupper($4), toupper($5), 4332, $7, $9}' >> gwas_pgc_bipolar_bmi_format.txt
head gwas_pgc_bipolar_bmi_format.txt

################################################################################
# Munge GWAS summary statistics for analysis in LDSC.

cd $path_gwas

# Accession of BMI GWAS summary statistics.
$path_ldsc/munge_sumstats.py \
--sumstats gwas_bmi_format.txt \
--out gwas_bmi \
--merge-alleles $path_alleles/w_hm3.snplist \

# BMI GWAS in Bipolar Disorder Controls in U.K. Biobank.
$path_ldsc/munge_sumstats.py \
--sumstats gwas_ukb_bipolar_control_bmi_format.txt \
--out gwas_ukb_bipolar_control_bmi \
--merge-alleles $path_alleles/w_hm3.snplist \

# BMI GWAS in Bipolar Disorder Cases in U.K. Biobank.
$path_ldsc/munge_sumstats.py \
--sumstats gwas_ukb_bipolar_case_bmi_format.txt \
--out gwas_ukb_bipolar_case_bmi \
--merge-alleles $path_alleles/w_hm3.snplist \

# Coombes BMI in Bipolar Disorder (PGC cohort) GWAS summary statistics.
#--signed-sumstats BETA,0 \ # I don't think this argument is consequential.
$path_ldsc/munge_sumstats.py \
--sumstats gwas_pgc_bipolar_bmi_format.txt \
--out gwas_pgc_bipolar_bmi \
--merge-alleles $path_alleles/w_hm3.snplist \

################################################################################
# Heritability in LDSC.


################################################################################
# Genetic correlation in LDSC.

if false; then
  $path_ldsc/ldsc.py \
  --rg gwas_bmi_raw.sumstats.gz,gwas_bipolar_bmi_coombes.sumstats.gz \
  --ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
  --w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
  --out correlation_bmi_control_raw_versus_bmi_bipolar_coombes.txt

  $path_ldsc/ldsc.py \
  --rg gwas_bmi_coombes.sumstats.gz,gwas_bipolar_bmi_coombes.sumstats.gz \
  --ref-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
  --w-ld-chr $path_disequilibrium/eur_w_ld_chr/ \
  --out correlation_bmi_control_coombes_versus_bmi_bipolar_coombes.txt
fi

################################################################################
# Deactivate Virtual Environment.

deactivate
which python2
