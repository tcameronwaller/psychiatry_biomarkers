#!/bin/bash

###########################################################################
# Organize paths.
# Read private, local file paths.
#echo "read private file path variables and organize paths..."
cd ~/paths
path_waller_tools=$(<"./waller_tools.txt")
path_environment_crossmap="${path_waller_tools}/python/environments/crossmap"
path_bcftools=$(<"./tools_bcftools.txt")
path_plink2=$(<"./tools_plink2.txt")
path_process=$(<"./process_psychiatric_metabolism.txt")
path_dock="${path_process}/dock"
path_parameters="${path_dock}/parameters"
path_translations_chromosomes_mayo="${path_parameters}/promiscuity/translations_chromosomes_mayo_bipolar.txt"
path_human_genome_sequence_compress="${path_dock}/access/human_genome_sequence/grch37/GRCh37.p13.genome.fa.gz"
path_human_genome_sequence="${path_dock}/access/human_genome_sequence/grch37/GRCh37.p13.genome.fa"
#host="ucsc"
host="ensembl"
if [[ "$host" == "ucsc" ]]; then
  path_assembly_translation_chain="${path_dock}/access/human_genome_assembly_chain/ucsc/hg38ToHg19.over.chain.gz"
elif [[ "$host" == "ensembl" ]]; then
  path_assembly_translation_chain="${path_dock}/access/human_genome_assembly_chain/ensembl/GRCh38_to_GRCh37.chain.gz"
else
  echo "invalid specification of host for genetic reference information"
fi
path_dbsnp="${path_dock}/access/dbsnp/grch37_format/GCF_000001405.25.gz"
path_mayo_bipolar_genotype_raw="${path_dock}/access/mayo_bipolar_genotype_raw"
path_mayo_bipolar_genotype_assembly="${path_dock}/access/mayo_bipolar_genotype_grch37_${host}"
path_mayo_bipolar_genotype_format="${path_dock}/access/mayo_bipolar_genotype_grch37_${host}_format"
path_genotype_snp_relevance_bim="${path_mayo_bipolar_genotype_format}/genotype_snp_relevance_bim"

# Scripts.
path_promiscuity_scripts="${path_process}/promiscuity/scripts"
path_script_submit_genotype_translate_assembly="${path_promiscuity_scripts}/utility/crossmap/1_submit_batch_directory_all_vcf_assembly_grch38_to_grch37.sh"
path_script_submit_genotype_format_annotation="${path_promiscuity_scripts}/utility/bcftools/1_submit_batch_directory_all_vcf_format_annotation.sh"
path_script_drive_extract_vcf_to_bim="${path_promiscuity_scripts}/utility/plink/drive_directory_all_extract_vcf_snps_to_plink_bim.sh"

###########################################################################
# Execute procedure.
###########################################################################

# Echo each command to console.
set -x

################################################################################
# Translate chromosome and base-pair position coordinates from human genome
# assembly GRCh38 to GRCh37.


path_test="${path_mayo_bipolar_genotype_raw}/MERGED.maf0.dosR20.3.noDups.chr21.dose.vcf.gz"
name_base_file="$(basename $path_test .vcf.gz)"
path_container="$(dirname $path_test)"

echo "base of file name: ${name_base_file}"
echo "path of directories to file: ${path_container}"

path_test="${path_dock}/test"
path_vcf_source="${path_mayo_bipolar_genotype_raw}/MERGED.maf0.dosR20.3.noDups.chr21.dose.vcf.gz"
path_vcf_product="${path_test}/MERGED.maf0.dosR20.3.noDups.chr21.dose.vcf.gz"

# Initialize directory.
rm -r $path_test
mkdir -p $path_test

threads=16
report="true"
path_script_translate_genome_assembly_vcf="${path_promiscuity_scripts}/utility/crossmap/translate_genome_assembly_vcf.sh"
/usr/bin/bash "${path_script_translate_genome_assembly_vcf}" \
$path_vcf_source \
$path_vcf_product \
$path_assembly_translation_chain \
$path_human_genome_sequence \
$threads \
$path_environment_crossmap \
$path_bcftools \
$report
