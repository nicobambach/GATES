#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# sourcing helper file to load run_cmd and log functions
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/scripts/helper.sh"

tumor_only_preprocess_workflow() {
    
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
align_tumor_only (){
    bwa mem \
    -R "@RG\tID:"$SAMPLE_NAME"\tPL:ILLUMINA\tSM:${SAMPLE_NAME}" \
    -t $THREADS \
    $REFERENCE \
    $TUMOR_FQ1 \
    $TUMOR_FQ2 \
    | samtools view -Sb - > ${MAPPED_READS_DIR}/${SAMPLE_NAME}.bam
}

run_cmd align_tumor_only

#### MARKING DUPLICATES AND SORTING ####
log "Marking duplicates and sorting..."
run_cmd gatk MarkDuplicatesSpark \
            -I ${MAPPED_READS_DIR}/${SAMPLE_NAME}.bam \
            -O ${MAPPED_READS_DIR}/${SAMPLE_NAME}_sorted_dedup.bam

# removing unsorted bam file to save space
rm ${MAPPED_READS_DIR}/${SAMPLE_NAME}.bam  

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

if [[ ! -f "${REFERENCE%.*}.dict" ]]; then
    run_cmd gatk CreateSequenceDictionary -R "$REFERENCE"
else
    log "FASTA sequence dictionary file already exists, skipping .dict creation for: $REFERENCE"
fi

log "Running BQSR..."

# running first pass (training)
run_cmd gatk BaseRecalibrator \
            -I ${MAPPED_READS_DIR}/${SAMPLE_NAME}_sorted_dedup.bam \
            -R $REFERENCE \
            --known-sites $BQSR_DB_DIR/Homo_sapiens_assembly38.dbsnp138.vcf \
            --known-sites $BQSR_DB_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz \
            --known-sites $BQSR_DB_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
            --interval-padding 100 \
            --intervals $INTERVAL_LIST \
            -O ${BQSR_OUT_DIR}/${SAMPLE_NAME}_recalibration_table.table

# running second step (recalibrating)
run_cmd gatk ApplyBQSR \
            -I ${MAPPED_READS_DIR}/${SAMPLE_NAME}_sorted_dedup.bam  \
            -R $REFERENCE \
            --interval-padding 100 \
            --intervals $INTERVAL_LIST \
            --bqsr-recal-file ${BQSR_OUT_DIR}/${SAMPLE_NAME}_recalibration_table.table \
            -O ${MAPPED_READS_DIR}/${SAMPLE_NAME}_sorted_dedup_recal.bam 

# running first pass again to be able to compare effect of bqsr
run_cmd gatk BaseRecalibrator \
            -I ${MAPPED_READS_DIR}/${SAMPLE_NAME}_sorted_dedup_recal.bam \
            -R $REFERENCE \
            --known-sites $BQSR_DB_DIR/Homo_sapiens_assembly38.dbsnp138.vcf \
            --known-sites $BQSR_DB_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz \
            --known-sites $BQSR_DB_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
            --interval-padding 100 \
            --intervals $INTERVAL_LIST \
            -O ${BQSR_OUT_DIR}/${SAMPLE_NAME}_post_recalibration_table.table

# analyzing effects of bqsr on base scoring
run_cmd gatk AnalyzeCovariates \
            -before ${BQSR_OUT_DIR}/${SAMPLE_NAME}_recalibration_table.table \
            -after ${BQSR_OUT_DIR}/${SAMPLE_NAME}_post_recalibration_table.table \
            -plots ${BQSR_OUT_DIR}/${SAMPLE_NAME}_recalibration_plots.pdf

# removing bqsr data and non-recalibrated bam file to save space 
run_cmd rm -rf "$BQSR_DB_DIR"
run_cmd rm -f ${MAPPED_READS_DIR}/${SAMPLE_NAME}_sorted_dedup.bam

}