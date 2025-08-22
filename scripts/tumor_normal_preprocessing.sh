#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# sourcing helper file to load run_cmd and log functions
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/scripts/helper.sh"


tumor_normal_preprocess_workflow(){

# making necessary directories
mkdir -p "$BQSR_DB_DIR" "$BQSR_OUT_DIR" "$MAPPED_READS_DIR"

#### ALIGNING TO REFERENCE ####
log "Indexing reference..."
if [[ ! -f "${REFERENCE}.amb" ]] || [[ ! -f "${REFERENCE}.ann" ]] || [[ ! -f "${REFERENCE}.bwt" ]] || [[ ! -f "${REFERENCE}.pac" ]] || [[ ! -f "${REFERENCE}.sa" ]]; then
    run_cmd bwa index "$REFERENCE"
else
    log "BWA index already exists, skipping indexing for: $REFERENCE"
fi

log "Aligning reads to reference..."
# aligning for tumor sample 
align_tumor(){
    bwa mem \
    -R "@RG\tID:${SAMPLE_NAME}\tPL:ILLUMINA\tSM:${SAMPLE_NAME}" \
    -t $THREADS \
    $REFERENCE \
    $TUMOR_FQ1 \
    $TUMOR_FQ2 \
    | samtools view -Sb - > ${MAPPED_READS_DIR}/${SAMPLE_NAME}_tumor.bam
}
run_cmd align_tumor

# aligning for normal sample
align_normal(){
    bwa mem \
    -R "@RG\tID:"$SAMPLE_NAME"\tPL:ILLUMINA\tSM:${SAMPLE_NAME}" \
    -t $THREADS \
    $REFERENCE \
    $NORMAL_FQ1 \
    $NORMAL_FQ2 \
    | samtools view -Sb - > ${MAPPED_READS_DIR}/${SAMPLE_NAME}_normal.bam
}
run_cmd align_normal

#### MARKING DUPLICATES AND SORTING ####
log "Marking duplicates and sorting..."

# marking duplicates and sorting on tumor 
run_cmd gatk MarkDuplicatesSpark \
        -I ${MAPPED_READS_DIR}/${SAMPLE_NAME}_tumor.bam \
        -O ${MAPPED_READS_DIR}/${SAMPLE_NAME}_tumor_sorted_dedup.bam

rm ${MAPPED_READS_DIR}/${SAMPLE_NAME}_tumor.bam 

# marking duplicagtes and sorting on normal 
run_cmd gatk MarkDuplicatesSpark \
        -I ${MAPPED_READS_DIR}/${SAMPLE_NAME}_normal.bam \
        -O ${MAPPED_READS_DIR}/${SAMPLE_NAME}_normal_sorted_dedup.bam

rm ${MAPPED_READS_DIR}/${SAMPLE_NAME}_normal.bam 

#### RUNNING BASE QUALITY SCORE RECALIBRATION ####
log "Downloading data for Base Quality Score Recalibration (BQSR)..."

#### RUNNING BASE QUALITY SCORE RECALIBRATION ####
log "Downloading data for Base Quality Score Recalibration (BQSR)..."

# getting 1000 genomes data
if [[ ! -f "$BQSR_DB_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" ]]; then
    run_cmd wget -P "$BQSR_DB_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
else
    log "File exists, skipping: Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
fi

if [[ ! -f "$BQSR_DB_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" ]]; then
    run_cmd wget -P "$BQSR_DB_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
else
    log "File exists, skipping: Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
fi

# getting known indels data
if [[ ! -f "$BQSR_DB_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz" ]]; then
    run_cmd wget -P "$BQSR_DB_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
else
    log "File exists, skipping: Homo_sapiens_assembly38.known_indels.vcf.gz"
fi

if [[ ! -f "$BQSR_DB_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi" ]]; then
    run_cmd wget -P "$BQSR_DB_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
else
    log "File exists, skipping: Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
fi

# getting dbsnp data
if [[ ! -f "$BQSR_DB_DIR/Homo_sapiens_assembly38.dbsnp138.vcf" ]]; then
    run_cmd wget -P "$BQSR_DB_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
else
    log "File exists, skipping: Homo_sapiens_assembly38.dbsnp138.vcf"
fi

if [[ ! -f "$BQSR_DB_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.idx" ]]; then
    run_cmd wget -P "$BQSR_DB_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
else
    log "File exists, skipping: Homo_sapiens_assembly38.dbsnp138.vcf.idx"
fi

# generating index file for reference fasta (needed for BQSR)
log "Generating reference FASTA index file"

if [[ ! -f "${REFERENCE}.fai" ]]; then
    run_cmd samtools faidx "$REFERENCE"
else
    log "FASTA index already exists, skipping indexing for: $REFERENCE"
fi

# generating dict file for reference fasta (needed for BQSR)
log "Generating reference FASTA sequence dictionary file"

if [[ ! -f "${REFERENCE}.dict" ]]; then
    run_cmd gatk CreateSequenceDictionary -R "$REFERENCE" -O "${REFERENCE}.dict"
else
    log "FASTA sequence dictionary file already exists, skipping .dict creation for: $REFERENCE"
fi

log "Running BQSR..."

# running bqsr on tumor (3 commands)
run_cmd gatk BaseRecalibrator \
    -I ${MAPPED_READS_DIR}/${SAMPLE_NAME}_tumor_sorted_dedup.bam \
    -R $REFERENCE \
    --known-sites $BQSR_DB_DIR/Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites $BQSR_DB_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites $BQSR_DB_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --interval-padding 100 \
    --intervals $INTERVAL_LIST \
    -O ${BQSR_OUT_DIR}/${SAMPLE_NAME}_tumor_recalibration_table.table

run_cmd gatk ApplyBQSR \
       -I ${MAPPED_READS_DIR}/${SAMPLE_NAME}_tumor_sorted_dedup.bam  \
       -R $REFERENCE \
       --interval-padding 100 \
       --intervals $INTERVAL_LIST \
       --bqsr-recal-file ${BQSR_OUT_DIR}/${SAMPLE_NAME}_tumor_recalibration_table.table \
       -O ${MAPPED_READS_DIR}/${SAMPLE_NAME}_tumor_sorted_dedup_recal.bam 

run_cmd gatk BaseRecalibrator \
    -I ${MAPPED_READS_DIR}/${SAMPLE_NAME}_tumor_sorted_dedup_recal.bam \
    -R $REFERENCE \
    --known-sites $BQSR_DB_DIR/Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites $BQSR_DB_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites $BQSR_DB_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --interval-padding 100 \
    --intervals $INTERVAL_LIST \
    -O ${BQSR_OUT_DIR}/${SAMPLE_NAME}_tumor_post_recalibration_table.table

# running bqsr on normal (3 commands)
run_cmd gatk BaseRecalibrator \
    -I ${MAPPED_READS_DIR}/${SAMPLE_NAME}_normal_sorted_dedup.bam \
    -R $REFERENCE \
    --known-sites $BQSR_DB_DIR/Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites $BQSR_DB_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites $BQSR_DB_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --interval-padding 100 \
    --intervals $INTERVAL_LIST \
    -O ${BQSR_OUT_DIR}/${SAMPLE_NAME}_normal_recalibration_table.table

run_cmd gatk ApplyBQSR \
       -I ${MAPPED_READS_DIR}/${SAMPLE_NAME}_normal_sorted_dedup.bam  \
       -R $REFERENCE \
       --interval-padding 100 \
       --intervals $INTERVAL_LIST \
       --bqsr-recal-file ${BQSR_OUT_DIR}/${SAMPLE_NAME}_normal_recalibration_table.table \
       -O ${MAPPED_READS_DIR}/${SAMPLE_NAME}_normal_sorted_dedup_recal.bam 

run_cmd gatk BaseRecalibrator \
    -I ${MAPPED_READS_DIR}/${SAMPLE_NAME}_normal_sorted_dedup_recal.bam \
    -R $REFERENCE \
    --known-sites $BQSR_DB_DIR/Homo_sapiens_assembly38.dbsnp138.vcf \
    --known-sites $BQSR_DB_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --known-sites $BQSR_DB_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --interval-padding 100 \
    --intervals $INTERVAL_LIST \
    -O ${BQSR_OUT_DIR}/${SAMPLE_NAME}_normal_post_recalibration_table.table

# analyzing effects of bqsr on tumor and normal samples
run_cmd gatk AnalyzeCovariates \
    -before ${BQSR_OUT_DIR}/${SAMPLE_NAME}_tumor_recalibration_table.table \
    -after ${BQSR_OUT_DIR}/${SAMPLE_NAME}_tumor_post_recalibration_table.table \
    -plots ${BQSR_OUT_DIR}/${SAMPLE_NAME}_tumor_recalibration_plots.pdf

run_cmd gatk AnalyzeCovariates \
    -before ${BQSR_OUT_DIR}/${SAMPLE_NAME}_normal_recalibration_table.table \
    -after ${BQSR_OUT_DIR}/${SAMPLE_NAME}_normal_post_recalibration_table.table \
    -plots ${BQSR_OUT_DIR}/${SAMPLE_NAME}_normal_recalibration_plots.pdf

# removing bqsr data and non-recalibrated bam files to save space 
run_cmd rm -rf "$BQSR_DB_DIR"
run_cmd rm -f ${MAPPED_READS_DIR}/${SAMPLE_NAME}_tumor_sorted_dedup.bam
run_cmd rm -f ${MAPPED_READS_DIR}/${SAMPLE_NAME}_normal_sorted_dedup.bam
}