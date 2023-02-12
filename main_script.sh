#!/bin/bash

### -------------- ###
### Initialization ###
### -------------- ###

# Loading modules
module load AWScli/aws-cli-2.1.0
module load sratoolkit.3.0.1
module load bwa/bwa-0.7.17
module load samtools/1.9
module load gatk/4.1.7.0

# Paths & Variables
reads_dir='/davidb/yatirsolan/NGS_variant_discovery/reads'
raw_reads1=${reads_dir}'/SRR5439568_1.fastq'
raw_reads2=${reads_dir}'/SRR5439568_2.fastq'
alignment_dir='/davidb/yatirsolan/NGS_variant_discovery/2nd_pipeline/alignment'
variants_dir='/davidb/yatirsolan/NGS_variant_discovery/2nd_pipeline/variants'
annotation_dir='/davidb/yatirsolan/NGS_variant_discovery/2nd_pipeline/annotation'
known_sites_dir='/davidb/yatirsolan/NGS_variant_discovery/2nd_pipeline/known_sites'
known_sites=${known_sites_dir}'/Homo_sapiens_assembly38.dbsnp138.vcf'
reference_genome_dir='/davidb/yatirsolan/NGS_variant_discovery/reference/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta'
reference_genome=${reference_genome_dir}'/Homo_sapiens_assembly38.fasta'
data_sources_dir='/davidb/yatirsolan/NGS_variant_discovery/2nd_pipeline/funcotator_dataSources.v1.6.20190124g'

# Creating directories 
mkdir /davidb/yatirsolan/NGS_variant_discovery
cd /davidb/yatirsolan/NGS_variant_discovery
mkdir ${reads_dir} ${alignment_dir} ${variants_dir} ${annotation_dir} ${known_sites_dir}

# Obtaining data
 # 1) Raw reads download from the NIH by the SRAtoolkit
fasterq-dump -O ${reads_dir} --split-files SRR5439568

 # 2) Reference genome donwload: 
aws s3 --no-sign-request --region eu-west-1 sync s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/ ./reference/Homo_sapiens/GATK/GRCh38/

 # 3) Downloading from the URL given below crucials files for BQSR stage.
    # https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?pli=1&prefix=&forceOnObjectsSortingFiltering=false
wget -P ${known_sites_dir} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget -P ${known_sites_dir} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

 # 4) Downloading a data sources for the functional annotation step 
gatk FuncotatorDataSourceDownloader -germline true -validate-integrity true -extract-after-download true 


### ------------------- ###
### Data pre-processing ###
### ------------------- ###

# Reads quality control
fastqc ${raw_reads1} ${raw_reads2}

# Reads mapping
read_group="@RG\tID:SRR5439568\tPL:ILLUMINA\tSM:KID"
bwa mem -R ${read_group} -t 4 ${reference_genome} ${raw_reads1} ${raw_reads2} | \ # mapping the reads to the reference genome.
samtools sort | \ # sorting the file in coordinance with the reference genome.
samtools view -b -h -o ${alignment_dir}/SRR5439568_srtd.bam # converting the final output to be in a BAM format.

# Flag duplicates reads 
gatk MarkDuplicates \
    -I ${alignment_dir}/SRR5439568_srtd.bam \
    -O ${alignment_dir}/SRR5439568_srtd_mrkdpl.bam \
    -M ${alignment_dir}/SRR5439568_srtd_mrkdpl.mtrx \
    -ASO coordinate # Assume sort order - coordinate, because the BAM file is already sorted in coordinate order.

# BQSR - Base quality scores adjacment 
 # 1) creating the model
gatk BaseRecalibrator \
    -R ${reference_genome} \
    -I ${alignment_dir}/SRR5439568_srtd_mrkdpl.bam \
    --known-sites ${known_sites} \
    -O ${alignment_dir}/recalibration.table
 # 2) applying the model
gatk ApplyBQSR \
    -R ${reference_genome} \
    -I ${alignment_dir}/SRR5439568_srtd_mrkdpl.bam \
    -bqsr ${alignment_dir}/recalibration.table \
    -O ${alignment_dir}/SRR5439568_srtd_mrkdpl_bqsr.bam

analaysis_ready_bam=${alignment_dir}/SRR5439568_srtd_mrkdpl_bqsr.bam

### ----------------- ###
### Variant Discovery ###
### ----------------- ###

# Variant calling
gatk HaplotypeCaller \
    -R ${reference_genome} \
    -I ${analaysis_ready_bam} \
    -O ${variants_dir}/SRR5439568_hptpcl.vcf

# Deviding the variants to SNOs and Indels
 # 1) SNPs
gatk SelectVariants \
    -V ${variants_dir}/SRR5439568_hptpcl.vcf \
    -select-type SNP \
    -O ${variants_dir}/SRR5439568_hptpcl_snps.vcf 
 # 2) Indels
gatk SelectVariants \
    -V ${variants_dir}/SRR5439568_hptpcl.vcf \
    -select-type INDEL \
    -O ${variants_dir}/SRR5439568_hptpcl_indels.vcf 

# Adding filtration tags
 # 1) SNPs
gatk VariantFiltration \
    -V ${variants_dir}/SRR5439568_hptpcl_snps.vcf \
    -O ${variants_dir}/SRR5439568_hptpcl_snps_tag.vcf \
    -filter-name 'QD_flt' -filter 'QD < 2.0' \
    -filter-name 'FS_flt' -filter 'FS > 60.0' \
    -filter-name 'SOR_flt' -filter 'SOR > 4.0' \
    -filter-name 'MQ_flt' -filter 'MQ < 40.0' \
    -filter-name 'MQRankSum_flt' -filter 'MQRankSum < -12.5' \
    -filter-name 'ReadPosRankSum_flt' -filter 'ReadPosRankSum < -8.0' \
    -genotype-filter-expression 'DP < 10' \
    -genotype-filter-name 'DP_flt' \
    -genotype-filter-expression 'GQ < 10' \
    -genotype-filter-name 'GQ_flt' 
 # 2) Indels
gatk VariantFiltration \
    -V ${variants_dir}/SRR5439568_hptpcl_indels.vcf \
    -O ${variants_dir}/SRR5439568_hptpcl_indels_tag.vcf \
    -filter-name 'QD_flt' -filter 'QD < 2.0' \
    -filter-name 'FS_flt' -filter 'FS > 200.0' \
    -filter-name 'SOR_flt' -filter 'SOR > 10.0' \
    -genotype-filter-expression 'DP < 10' \
    -genotype-filter-name 'DP_flt' \
    -genotype-filter-expression 'GQ < 10' \
    -genotype-filter-name 'GQ_flt' 

# Excluding the filtered variants
 # 1.1) SNPs - Exclude variants that failed to pass site-level filtration.
gatk SelectVariants \
    -V ${variants_dir}/SRR5439568_hptpcl_snps_tag.vcf \
    -O ${variants_dir}/SRR5439568_hptpcl_snps_fltrd.vcf \
    --exclude-filtered true 
 # 1.2) SNPs - Exclude variants that failed to pass sample-level filtration.
echo "$(grep -v -E 'DP_flt|GQ_flt' ${variants_dir}/SRR5439568_hptpcl_snps_fltrd.vcf)" > \
      ${variants_dir}/SRR5439568_hptpcl_snps_fltrd.vcf

# --------------- for the indels, the same script except file names --------------- #

 # 2.1) Indels - Exclude variants that failed to pass site-level filtration.
gatk SelectVariants \
    -V ${variants_dir}/SRR5439568_hptpcl_indels_tag.vcf \
    -O ${variants_dir}/SRR5439568_hptpcl_indels_fltrd.vcf \
    --exclude-filtered true
 # 2.2) Indels - Exclude variants that failed to pass sample-level filtration.
echo "$(grep -v -E 'DP_flt|GQ_flt' ${variants_dir}/SRR5439568_hptpcl_indels_fltrd.vcf)" > \
      ${variants_dir}/SRR5439568_hptpcl_indels_fltrd.vcf

annotation_ready_SNPs=${variants_dir}/SRR5439568_hptpcl_snps_fltrd.vcf
annotation_ready_indels=${variants_dir}/SRR5439568_hptpcl_indels_fltrd.vcf

### ------------------ ###
### Variant annotation ### 
### ------------------ ###

# Adding a functional annotation to the INFO column of the VCF.
 # 1) SNPs
gatk Funcotator \
    -R ${reference_genome} \
    --ref-version hg38 \
    --data-sources-path ${data_sources_dir} \
    -V ${annotation_ready_SNPs} \
    -O ${annotation_dir}/SRR5439568_hptpcl_snps_fltrd_anttd.vcf \
    --output-file-format VCF
# --------------- for the indels, the same script except file names --------------- #
 # 2) Indels
gatk Funcotator \
    -R ${reference_genome} \
    --ref-version hg38 \
    --data-sources-path ${data_sources_dir} \
    -V ${annotation_ready_indels} \
    -O ${annotation_dir}/SRR5439568_hptpcl_indels_fltrd_anttd.vcf \
    --output-file-format VCF

# Creating a table with variants and their annotations
 # 1) SNPs 
  # 1.1) Extracting the column's names from the header section of the VCF into an a new empty file. 
funcotation_header=$(grep '^##' ${annotation_dir}/SRR5439568_hptpcl_snps_fltrd_anttd.vcf | \
                   grep 'FUNCOTATION' | \
                   sed 's/ //g;s/>//g;s/"//g') # remove spaces, double quotes, and right arrows from the string.
echo ${funcotation_header##*\Funcotationfieldsare:} | \
     sed 's/|/\t/g' > \ # replace pipes '|' with tabs '\t'.
     ${annotation_dir}/SRR5439568_annotated_variants_snps.tsv

  # 1.2) Extracting the INFO column from the VCF, with the FUNCOTATION data only, and dropping it to a temporary file.
gatk VariantsToTable \
    -V ${annotation_dir}/SRR5439568_hptpcl_snps_fltrd_anttd.vcf \
    -O ${annotation_dir}/SRR5439568_annotated_variants_snps_tmp.tsv \
    -F FUNCOTATION 
    
  # 1.3) Uniting the column's name and the data itself to form a table (tsv file).
grep -v 'FUNCOTATION' ${annotation_dir}/SRR5439568_annotated_variants_snps_tmp.tsv | \
sed 's/[][]//g;s/|/\t/g' \ # removing square brackets from the string, and replace pipes '|' with tabs '\t'.
>> ${annotation_dir}/SRR5439568_annotated_variants_snps.tsv

rm ${annotation_dir}/SRR5439568_annotated_variants_snps_tmp.tsv

# --------------- for the indels, the same script except file names --------------- #

 # 2) Indels
  # 2.1) Extracting the column's names from the header section of the VCF into an a new empty file. 
funcotation_header=$(grep '^##' ${annotation_dir}/SRR5439568_hptpcl_indels_fltrd_anttd.vcf | \
                   grep 'FUNCOTATION' | \
                   sed 's/ //g;s/>//g;s/"//g') # remove spaces, double quotes, and right arrows from the string.
echo ${funcotation_header##*\Funcotationfieldsare:} | \
     sed 's/|/\t/g' > \ # replace pipes '|' with tabs '\t'.
     ${annotation_dir}/SRR5439568_annotated_variants_indels.tsv # adding the columns to the head of an empty file

  # 2.2) Extracting the INFO column from the VCF, with the FUNCOTATION data only, and dropping it to a temporary file.
gatk VariantsToTable \
    -V ${annotation_dir}/SRR5439568_hptpcl_indels_fltrd_anttd.vcf \
    -O ${annotation_dir}/SRR5439568_annotated_variants_indels_tmp.tsv \
    -F FUNCOTATION 

  # 2.3) Uniting the column's name and the data itself to form a table (tsv file).
grep -v "FUNCOTATION" ${annotation_dir}/SRR5439568_annotated_variants_indels_tmp.tsv | \
sed 's/[][]//g;s/|/\t/g' \ # removing square brackets from the string, and replace pipes '|' with tabs '\t'.
>> ${annotation_dir}/SRR5439568_annotated_variants_indels.tsv

rm ${annotation_dir}/SRR5439568_annotated_variants_indels_tmp.tsv # adding the columns to the head of an empty file

###
### --------- the end --------- ###
###