#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# getting root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# sourcing helper file to load run_cmd and log functions
source "$SCRIPT_DIR/scripts/helper.sh"

germline_variant_calling_workflow(){

# setting directories
ANALYSIS_DIR="./germline"
HAPLOTYPECALLER_OUTPUT_DIR="$ANALYSIS_DIR/variants"
REALIGNED_READS_DIR="$ANALYSIS_DIR/realigned_reads"

# making necessary directories
mkdir -p "$HAPLOTYPECALLER_OUTPUT_DIR" "$REALIGNED_READS_DIR"

# getting sample name from BAM header
SAMPLE_NAME=$(samtools samples $TUMOR_BAM | cut -f1)

log "Beginning germline variant calling for $SAMPLE_NAME"

#### RUNNING HAPLOTYPECALLER ####
log "Running HaplotypeCaller..."
run_cmd gatk HaplotypeCaller \
    -R $REFERENCE \
    -I $TUMOR_BAM \
    -L $INTERVAL_LIST \
    --native-pair-hmm-threads $THREADS \
    -O ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_variants.vcf.gz \
    --bam-output ${REALIGNED_READS_DIR}/${SAMPLE_NAME}_realigned.bam

#### SPLITTING SNPS AND INDELS FOR HARD-FILTERING ####
log "Splitting SNP and INDELs for hard-filtering..."

# selecting SNPs
run_cmd gatk SelectVariants \
    -R $REFERENCE \
    -V ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_variants.vcf.gz \
    -select-type SNP \
    -O ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_snps.vcf.gz

# selecting INDELs
run_cmd gatk SelectVariants \
    -R $REFERENCE \
    -V ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_variants.vcf.gz \
    -select-type INDEL \
    -O ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_indels.vcf.gz

# removing non-filtered combined VCF
rm -f ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_variants.vcf.gz
rm -f ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_variants.vcf.gz.tbi

#### APPLYING HARD FILTERING TO SNP and INDEL VCF FILES ####
log "Applying standard hard-filters for SNPs and INDELs..."

# filtration for SNPs
run_cmd gatk VariantFiltration \
    -R $REFERENCE \
    -V ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRS-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "RPRS-8" \
    -O ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_all_snps.vcf.gz

# removing non-filtered SNP VCFs
rm -f ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_snps.vcf.gz
rm -f ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_snps.vcf.gz.tbi

# filtration for INDELs
run_cmd gatk VariantFiltration \
    -R $REFERENCE \
    -V ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_indels.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "SOR > 10.0" --filter-name "SOR10" \
    -filter "ReadPosRankSum < -20.0" --filter-name "RPRS-20" \
    -O ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_all_indels.vcf.gz

# removing non-filtered INDEL VCFs
rm -f ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_indels.vcf.gz
rm -f ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_indels.vcf.gz.tbi

#### COMBINING SNP AND INDEL VCF FILES ####
log "Combining filtered SNP and INDEL VCF files..."

run_cmd gatk MergeVcfs \
    -I ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_all_snps.vcf.gz \
    -I ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_all_indels.vcf.gz \
    -O ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_all_germline_variants.vcf.gz

# removing non-combined VCF files
rm -f ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_all_snps.vcf.gz
rm -f ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_all_snps.vcf.gz.tbi
rm -f ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_all_indels.vcf.gz
rm -f ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_all_indels.vcf.gz.tbi

#### SELECTING ONLY PASSED VARIANTS ####
log "Selecting variants that passed filtering..."

run_cmd gatk SelectVariants \
    -R $REFERENCE \
    -V ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_all_germline_variants.vcf.gz \
    -O ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_passed_germline_variants.vcf.gz \
    --exclude-filtered true

#### FINAL STATISTICS ####
ALL_VARIANTS=$(bcftools view -H ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_all_germline_variants.vcf.gz | wc -l | tr -d ' ')
PASSED_VARIANTS=$(bcftools view -H ${HAPLOTYPECALLER_OUTPUT_DIR}/${SAMPLE_NAME}_passed_germline_variants.vcf.gz | wc -l | tr -d ' ')
FILTERED_OUT=$((ALL_VARIANTS - PASSED_VARIANTS))
PASS_PERCENT=$(( (PASSED_VARIANTS * 100) / ALL_VARIANTS ))

log "Finished germline variant calling for $SAMPLE_NAME"
log "Total variants found: ${ALL_VARIANTS}"
log "Variants filtered out: ${FILTERED_OUT}"
log "Variants that passed filtering: ${PASSED_VARIANTS} (${PASS_PERCENT}%)"
log "Variant output files:"
log "  - ${SAMPLE_NAME}_all_germline_variants.vcf.gz (all variants)"
log "  - ${SAMPLE_NAME}_passed_germline_variants.vcf.gz (variants that passed filtering)"

}