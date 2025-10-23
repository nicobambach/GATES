#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# getting root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# sourcing helper file to load run_cmd and log functions
source "$SCRIPT_DIR/scripts/helper.sh"

tumor_normal_variant_calling_workflow() {

# setting directories
SUPPORTING_FILES_DIR="./supporting_files"
MUTECT2_SUPPORTING_FILES_DIR="$SUPPORTING_FILES_DIR/calling_resources"
ANALYSIS_DIR="./tumor_normal_somatic"
MUTECT2_OUTPUT_DIR="$ANALYSIS_DIR/variants"
MUTECT2_FILTERING_DIR="$ANALYSIS_DIR/filtering_data"

# getting sample names from BAM headers
TUMOR_NAME=$(samtools samples $TUMOR_BAM | cut -f1)
NORMAL_NAME=$(samtools samples $NORMAL_BAM | cut -f1)

log "Beginning tumor-normal somatic variant calling for $TUMOR_NAME"
log "Tumor sample: $TUMOR_NAME"
log "Normal sample: $NORMAL_NAME"

# making necessary directories
mkdir -p "$MUTECT2_SUPPORTING_FILES_DIR" "$MUTECT2_OUTPUT_DIR" "$MUTECT2_FILTERING_DIR"

#### DOWNLOADING NECESSARY MUTECT2 PON AND GENOMEAD FILES ####
log "Downloading necessary supporting files for Mutect2..."

# getting public panel of normals (PON) from GATK resources
if [[ ! -f "$MUTECT2_SUPPORTING_FILES_DIR"/1000g_pon.hg38.vcf.gz ]]; then
    run_cmd wget -P "$MUTECT2_SUPPORTING_FILES_DIR" https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
else
    log "File exists, skipping download for: 1000g_pon.hg38.vcf.gz"
fi

if [[ ! -f "$MUTECT2_SUPPORTING_FILES_DIR"/1000g_pon.hg38.vcf.gz.tbi ]]; then
    run_cmd wget -P "$MUTECT2_SUPPORTING_FILES_DIR" https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi
else
    log "File exists, skipping download for: 1000g_pon.hg38.vcf.gz.tbi"
fi

# getting annotated germline variants from GATK resoucres
if [[ ! -f "$MUTECT2_SUPPORTING_FILES_DIR"/af-only-gnomad.hg38.vcf.gz ]]; then 
    run_cmd wget -P "$MUTECT2_SUPPORTING_FILES_DIR" https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
else 
    log "File exists, skipping download for: af-only-gnomad.hg38.vcf.gz"
fi

if [[ ! -f "$MUTECT2_SUPPORTING_FILES_DIR"/af-only-gnomad.hg38.vcf.gz.tbi ]]; then 
    run_cmd wget -P "$MUTECT2_SUPPORTING_FILES_DIR" https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
else 
    log "File exists, skipping download for: af-only-gnomad.hg38.vcf.gz.tbi"
fi 

#### RUNNING MUTECT2 IN TUMOR-NORMAL MODE ####
log "Running Mutect2 in tumor-normal mode..."

run_cmd gatk Mutect2 \
    -I $TUMOR_BAM \
    -tumor $TUMOR_NAME \
    -I $NORMAL_BAM \
    -normal $NORMAL_NAME \
    -R $REFERENCE \
    --germline-resource $MUTECT2_SUPPORTING_FILES_DIR/af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals $MUTECT2_SUPPORTING_FILES_DIR/1000g_pon.hg38.vcf.gz \
    --intervals $INTERVAL_LIST \
    --interval-padding 0 \
    --f1r2-tar-gz ${MUTECT2_FILTERING_DIR}/${TUMOR_NAME}_f1r2.tar.gz \
    --native-pair-hmm-threads $THREADS \
    -O ${MUTECT2_OUTPUT_DIR}/${TUMOR_NAME}_variants.vcf.gz

#### ESTIMATING CROSS-SAMPLE CONTAMINATION ####
log "Downloading common SNP reference VCF file..."

if [[ ! -f "$MUTECT2_SUPPORTING_FILES_DIR"/small_exac_common_3.hg38.vcf.gz ]]; then
    run_cmd wget -P "$MUTECT2_SUPPORTING_FILES_DIR" https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz
else
    log "File exists, skipping: small_exac_common_3.hg38.vcf.gz"
fi

if [[ ! -f "$MUTECT2_SUPPORTING_FILES_DIR"/small_exac_common_3.hg38.vcf.gz.tbi ]]; then
    run_cmd wget -P "$MUTECT2_SUPPORTING_FILES_DIR" https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi
else
    log "File exists, skipping: small_exac_common_3.hg38.vcf.gz.tbi"
fi

log "Estimating cross-sample contamination..."

# getting pileups at common SNP sites for tumor
run_cmd gatk GetPileupSummaries \
    -I $TUMOR_BAM \
    -V $MUTECT2_SUPPORTING_FILES_DIR/small_exac_common_3.hg38.vcf.gz \
    -L $MUTECT2_SUPPORTING_FILES_DIR/small_exac_common_3.hg38.vcf.gz \
    -L $INTERVAL_LIST \
    -O ${MUTECT2_FILTERING_DIR}/${TUMOR_NAME}_pileups_summary.table

# getting pileups at common SNP sites for normal
run_cmd gatk GetPileupSummaries \
    -I $NORMAL_BAM \
    -V $MUTECT2_SUPPORTING_FILES_DIR/small_exac_common_3.hg38.vcf.gz \
    -L $MUTECT2_SUPPORTING_FILES_DIR/small_exac_common_3.hg38.vcf.gz \
    -L $INTERVAL_LIST \
    -O ${MUTECT2_FILTERING_DIR}/${NORMAL_NAME}_pileups_summary.table

# running CalculateContamination with tumor and normal
run_cmd gatk CalculateContamination \
    -I ${MUTECT2_FILTERING_DIR}/${TUMOR_NAME}_pileups_summary.table \
    -matched ${MUTECT2_FILTERING_DIR}/${NORMAL_NAME}_pileups_summary.table \
    -tumor-segmentation ${MUTECT2_FILTERING_DIR}/${TUMOR_NAME}_segmentation_table.table \
    -O ${MUTECT2_FILTERING_DIR}/${TUMOR_NAME}_contamination_table.table

#### ESTIMATING READ ORIENTATION BIAS ####
log "Estimating read orientation bias..."

run_cmd gatk LearnReadOrientationModel \
    -I ${MUTECT2_FILTERING_DIR}/${TUMOR_NAME}_f1r2.tar.gz \
    -O ${MUTECT2_FILTERING_DIR}/${TUMOR_NAME}_read_orientation_model.tar.gz

#### FILTERING VARIANT CALLS ####
log "Filtering called variants..."

run_cmd gatk FilterMutectCalls \
    -R $REFERENCE \
    -V ${MUTECT2_OUTPUT_DIR}/${TUMOR_NAME}_variants.vcf.gz  \
    --contamination-table ${MUTECT2_FILTERING_DIR}/${TUMOR_NAME}_contamination_table.table \
    --tumor-segmentation ${MUTECT2_FILTERING_DIR}/${TUMOR_NAME}_segmentation_table.table \
    --orientation-bias-artifact-priors ${MUTECT2_FILTERING_DIR}/${TUMOR_NAME}_read_orientation_model.tar.gz \
    --filtering-stats ${MUTECT2_OUTPUT_DIR}/${TUMOR_NAME}_filtering_stats.txt \
    -O ${MUTECT2_OUTPUT_DIR}/${TUMOR_NAME}_all_somatic_variants.vcf.gz

# removing non-filtered VCF files
rm -f ${MUTECT2_OUTPUT_DIR}/${TUMOR_NAME}_variants.vcf.gz  
rm -f ${MUTECT2_OUTPUT_DIR}/${TUMOR_NAME}_variants.vcf.gz.tbi
rm -f ${MUTECT2_OUTPUT_DIR}/${TUMOR_NAME}_variants.vcf.gz.stats

log "Selecting only variants that passed filtering..."
run_cmd gatk SelectVariants \
    -V ${MUTECT2_OUTPUT_DIR}/${TUMOR_NAME}_all_somatic_variants.vcf.gz \
    --exclude-filtered true \
    -O ${MUTECT2_OUTPUT_DIR}/${TUMOR_NAME}_passed_somatic_variants.vcf.gz

#### FINAL STATISTICS ####
ALL_VARIANTS=$(bcftools view -H ${MUTECT2_OUTPUT_DIR}/${TUMOR_NAME}_all_somatic_variants.vcf.gz | wc -l)
PASSED_VARIANTS=$(bcftools view -H ${MUTECT2_OUTPUT_DIR}/${TUMOR_NAME}_passed_somatic_variants.vcf.gz | wc -l)
FILTERED_OUT=$((ALL_VARIANTS - PASSED_VARIANTS))
PASS_PERCENT=$(( (PASSED_VARIANTS * 100) / ALL_VARIANTS ))

log "Finished tumor-normal somatic variant calling for $TUMOR_NAME"
log "Total variants called: ${ALL_VARIANTS}"
log "Variants filtered out: ${FILTERED_OUT}"
log "Variants that passed filtering: ${PASSED_VARIANTS} (${PASS_PERCENT}%)"
log "Variant output files:"
log "  - ${TUMOR_NAME}_all_somatic_variants.vcf.gz (all variants)"
log "  - ${TUMOR_NAME}_passed_somatic_variants.vcf.gz (variants that passed filtering)"

}