#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# sourcing helper file to load run_cmd and log functions
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/scripts/helper.sh"

preprocess_workflow() {

# setting directories
SUPPORTING_FILES_DIR="./supporting_files"
PREPROCESSING_DIR="./preprocessing"
BQSR_DB_DIR="$SUPPORTING_FILES_DIR/preprocessing_resources"
BQSR_OUT_DIR="$PREPROCESSING_DIR/bqsr_output"
MAPPED_READS_DIR="$PREPROCESSING_DIR/mapped_reads"
QC_DIR="$PREPROCESSING_DIR/qc"
TRIMMED_DIR="$PREPROCESSING_DIR/trimmed_reads" 

log "Starting preprocessing for $SAMPLE_NAME"
    
# making necessary directories
mkdir -p "$BQSR_DB_DIR" "$BQSR_OUT_DIR" "$MAPPED_READS_DIR" "$QC_DIR" "$TRIMMED_DIR"

#### FASTP QUALITY CONTROL AND TRIMMING ####
log "Running QC and adapter trimming..."

run_cmd fastp \
    -i "$FQ1" \
    -I "$FQ2" \
    -o "${TRIMMED_DIR}/${SAMPLE_NAME}_R1_trimmed.fastq.gz" \
    -O "${TRIMMED_DIR}/${SAMPLE_NAME}_R2_trimmed.fastq.gz" \
    --thread $THREADS \
    --detect_adapter_for_pe \
    --html "${QC_DIR}/${SAMPLE_NAME}_fastp.html" \
    --json "${QC_DIR}/${SAMPLE_NAME}_fastp.json"

FQ1="${TRIMMED_DIR}/${SAMPLE_NAME}_R1_trimmed.fastq.gz"
FQ2="${TRIMMED_DIR}/${SAMPLE_NAME}_R2_trimmed.fastq.gz" 

#### ALIGNING TO REFERENCE ####
log "Indexing reference..."
if [[ ! -f "${REFERENCE}.amb" ]] || [[ ! -f "${REFERENCE}.ann" ]] || [[ ! -f "${REFERENCE}.bwt" ]] || [[ ! -f "${REFERENCE}.pac" ]] || [[ ! -f "${REFERENCE}.sa" ]]; then
    run_cmd bwa index "$REFERENCE"
else
    log "BWA index already exists, skipping indexing for: $REFERENCE"
fi

log "Aligning reads to reference..."
align_sample (){
    bwa mem \
    -R "@RG\tID:"$SAMPLE_NAME"\tPL:ILLUMINA\tSM:${SAMPLE_NAME}" \
    -t $THREADS \
    $REFERENCE \
    $FQ1 \
    $FQ2 \
    | samtools view -Sb - > ${MAPPED_READS_DIR}/${SAMPLE_NAME}.bam
}

run_cmd align_sample

#### MARKING DUPLICATES AND SORTING ####
log "Marking duplicates and sorting..."
run_cmd gatk MarkDuplicatesSpark \
            -I ${MAPPED_READS_DIR}/${SAMPLE_NAME}.bam \
            -O ${MAPPED_READS_DIR}/${SAMPLE_NAME}_mrkdp.bam

# removing unsorted bam file to save space
rm ${MAPPED_READS_DIR}/${SAMPLE_NAME}.bam  

#### RUNNING BASE QUALITY SCORE RECALIBRATION ####
log "Downloading data for Base Quality Score Recalibration (BQSR)..."

# getting 1000 genomes data
if [[ ! -f "$BQSR_DB_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" ]]; then
    run_cmd wget -P "$BQSR_DB_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
else
    log "File exists, skipping download for: Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
fi

if [[ ! -f "$BQSR_DB_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" ]]; then
    run_cmd wget -P "$BQSR_DB_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
else
    log "File exists, skipping download for: Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
fi

# getting known indels data
if [[ ! -f "$BQSR_DB_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz" ]]; then
    run_cmd wget -P "$BQSR_DB_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
else
    log "File exists, skipping download for: Homo_sapiens_assembly38.known_indels.vcf.gz"
fi

if [[ ! -f "$BQSR_DB_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi" ]]; then
    run_cmd wget -P "$BQSR_DB_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
else
    log "File exists, skipping download for: Homo_sapiens_assembly38.known_indels.vcf.gz.tbi"
fi

# getting dbsnp data
if [[ ! -f "$BQSR_DB_DIR/Homo_sapiens_assembly38.dbsnp138.vcf" ]]; then
    run_cmd wget -P "$BQSR_DB_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
else
    log "File exists, skipping download for: Homo_sapiens_assembly38.dbsnp138.vcf"
fi

if [[ ! -f "$BQSR_DB_DIR/Homo_sapiens_assembly38.dbsnp138.vcf.idx" ]]; then
    run_cmd wget -P "$BQSR_DB_DIR" https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
else
    log "File exists, skipping download for: Homo_sapiens_assembly38.dbsnp138.vcf.idx"
fi

# generating index file for reference fasta (needed for BQSR)
log "Generating reference FASTA index file..."

if [[ ! -f "${REFERENCE}.fai" ]]; then
    run_cmd samtools faidx "$REFERENCE"
else
    log "FASTA index already exists, skipping indexing for: $REFERENCE"
fi

# generating dict file for reference fasta (needed for BQSR)
log "Generating reference FASTA sequence dictionary file..."

if [[ ! -f "${REFERENCE%.*}.dict" ]]; then
    run_cmd gatk CreateSequenceDictionary -R "$REFERENCE"
else
    log "FASTA sequence dictionary file already exists, skipping .dict creation for: $REFERENCE"
fi

log "Running BQSR..."

# running first pass (training)
run_cmd gatk BaseRecalibrator \
            -I ${MAPPED_READS_DIR}/${SAMPLE_NAME}_mrkdp.bam \
            -R $REFERENCE \
            --known-sites $BQSR_DB_DIR/Homo_sapiens_assembly38.dbsnp138.vcf \
            --known-sites $BQSR_DB_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz \
            --known-sites $BQSR_DB_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
            --interval-padding 0 \
            --intervals $INTERVAL_LIST \
            -O ${BQSR_OUT_DIR}/${SAMPLE_NAME}_recal_tbl.table

# running second step (recalibrating)
run_cmd gatk ApplyBQSR \
            -I ${MAPPED_READS_DIR}/${SAMPLE_NAME}_mrkdp.bam  \
            -R $REFERENCE \
            --interval-padding 0 \
            --intervals $INTERVAL_LIST \
            --bqsr-recal-file ${BQSR_OUT_DIR}/${SAMPLE_NAME}_recal_tbl.table \
            -O ${MAPPED_READS_DIR}/${SAMPLE_NAME}_recal.bam 

# removing non-recalibrated bam file to save space 
run_cmd rm -f ${MAPPED_READS_DIR}/${SAMPLE_NAME}_mrkdp.bam
run_cmd rm -f ${MAPPED_READS_DIR}/${SAMPLE_NAME}_mrkdp.bam.bai
run_cmd rm -f ${MAPPED_READS_DIR}/${SAMPLE_NAME}_mrkdp.bam.sbi

#### COMPILING QC METRICS ####
log "Compiling additional quality control metrics..."

run_cmd samtools stats ${MAPPED_READS_DIR}/${SAMPLE_NAME}_recal.bam > ${QC_DIR}/${SAMPLE_NAME}_alignment_stats.txt
run_cmd mosdepth \
    --threads $THREADS \
    --by $INTERVAL_LIST \
    --mapq 20 \
    --no-per-base \
    --fast-mode \
    ${QC_DIR}/${SAMPLE_NAME} \
    ${MAPPED_READS_DIR}/${SAMPLE_NAME}_recal.bam

run_cmd multiqc ${QC_DIR} -f --outdir ${QC_DIR} \
    --filename ${SAMPLE_NAME}_multiqc_report.html \
    --title "QC Report for: ${SAMPLE_NAME}" \
    --force

log "Preprocessing complete for $SAMPLE_NAME" 
log "Final BAM: ${MAPPED_READS_DIR}/${SAMPLE_NAME}_recal.bam"
log "QC reports: ${QC_DIR}"

}