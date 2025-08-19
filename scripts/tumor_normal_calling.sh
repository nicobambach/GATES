#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# getting root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# sourcing helper file to load run_cmd and log functions
source "$SCRIPT_DIR/scripts/helper.sh"

tumor_normal_variant_calling_workflow() {

# making necessary directories
mkdir -p "$MUTECT2_SUPPORTING_FILES_DIR" "$MUTECT2_OUTPUT_DIR" "$MUTECT2_FILTERING_DIR"

#### DOWNLOADING NECESSARY MUTECT2 PON AND GENOMEAD FILES ####
log "Downloading necessary supporting files for Mutect2..."

# getting public panel of normals (PON) from GATK resources
if [[ ! -f "$MUTECT2_SUPPORTING_FILES_DIR"/1000g_pon.hg38.vcf.gz ]]; then
    run_cmd wget -P "$MUTECT2_SUPPORTING_FILES_DIR" https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
else
    log "File exists, skipping: 1000g_pon.hg38.vcf.gz"
fi

if [[ ! -f "$MUTECT2_SUPPORTING_FILES_DIR"/1000g_pon.hg38.vcf.gz.tbi ]]; then
    run_cmd wget -P "$MUTECT2_SUPPORTING_FILES_DIR" https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi
else
    log "File exists, skipping: 1000g_pon.hg38.vcf.gz.tbi"
fi

# getting annotated germline variants from GATK resoucres
if [[ ! -f "$MUTECT2_SUPPORTING_FILES_DIR"/af-only-gnomad.hg38.vcf.gz ]]; then 
    run_cmd wget -P "$MUTECT2_SUPPORTING_FILES_DIR" https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
else 
    log "File exists, skipping: af-only-gnomad.hg38.vcf.gz"
fi

if [[ ! -f "$MUTECT2_SUPPORTING_FILES_DIR"/af-only-gnomad.hg38.vcf.gz.tbi ]]; then 
    run_cmd wget -P "$MUTECT2_SUPPORTING_FILES_DIR" https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
else 
    log "File exists, skipping: af-only-gnomad.hg38.vcf.gz.tbi"
fi 

#### RUNNING MUTECT2 IN TUMOR-NORMAL MODE ####
log "Running Mutect2 in tumor-normal mode..."

run_cmd gatk Mutect2 \
    -I $TUMOR_BAM \
    -tumor TUMOR \
    -I $NORMAL_BAM \
    -normal NORMAL \
    -R $REFERENCE \
    --germline-resource $MUTECT2_SUPPORTING_FILES_DIR/af-only-gnomad.hg38.vcf.gz \
    --panel-of-normals $MUTECT2_SUPPORTING_FILES_DIR/1000g_pon.hg38.vcf.gz \
    --intervals $INTERVAL_LIST \
    --interval-padding 100 \
    --f1r2-tar-gz ${MUTECT2_OUTPUT_DIR}/${SAMPLE_NAME}_f1r2.tar.gz \
    -O ${MUTECT2_OUTPUT_DIR}/${SAMPLE_NAME}_variants.vcf.gz

 #### ESTIMATING CROSS-SAMPLE CONTAMINATION ####
log "Estimating cross-sample contamination..."

# running GetPileupSummaries for tumor
run_cmd gatk GetPileupSummaries \
    -I $TUMOR_BAM \
    -V $MUTECT2_SUPPORTING_FILES_DIR/af-only-gnomad.hg38.vcf.gz \
    --intervals $INTERVAL_LIST \
    --interval-padding 100 \
    -O ${MUTECT2_FILTERING_DIR}/${SAMPLE_NAME}_TUMOR_pileups_summary.table

# running GetPileupSummaries for normal
run_cmd gatk GetPileupSummaries \
    -I $NORMAL_BAM \
    -V $MUTECT2_SUPPORTING_FILES_DIR/af-only-gnomad.hg38.vcf.gz \
    --intervals $INTERVAL_LIST \
    --interval-padding 100 \
    -O ${MUTECT2_FILTERING_DIR}/${SAMPLE_NAME}_NORMAL_pileups_summary.table

# running CalculateContamination with 
run_cmd gatk CalculateContamination \
    -I ${MUTECT2_FILTERING_DIR}/${SAMPLE_NAME}_TUMOR_pileups_summary.table \
    -matched ${MUTECT2_FILTERING_DIR}/${SAMPLE_NAME}_NORMAL_pileups_summary.table \
    -tumor-segmentation ${MUTECT2_FILTERING_DIR}/${SAMPLE_NAME}_tumor_segmentation_table.table \
    -O ${MUTECT2_FILTERING_DIR}/${SAMPLE_NAME}_contamination_table.table

#### ESTIMATING READ ORIENTATION BIAS ####
log "Estimating read orientation bias..."

run_cmd gatk LearnReadOrientationModel \
    -I ${MUTECT2_OUTPUT_DIR}/${SAMPLE_NAME}_f1r2.tar.gz \
    -O ${MUTECT2_FILTERING_DIR}/${SAMPLE_NAME}_read_orientation_model.tar.gz

#### FILTERING VARIANT CALLS ####
log "Filtering called variants..."

run_cmd gatk FilterMutectCalls \
    -R $REFERENCE \
    -V ${MUTECT2_OUTPUT_DIR}/${SAMPLE_NAME}_variants.vcf.gz  \
    --contamination-table ${MUTECT2_FILTERING_DIR}/${SAMPLE_NAME}_contamination_table.table \
    --tumor-segmentation ${MUTECT2_FILTERING_DIR}/${SAMPLE_NAME}_tumor_segmentation_table.table \
    --orientation-bias-artifact-priors ${MUTECT2_FILTERING_DIR}/${SAMPLE_NAME}_read_orientation_model.tar.gz \
    --filtering-stats ${MUTECT2_OUTPUT_DIR}/${SAMPLE_NAME}_filtering_stats.txt \
    -O ${MUTECT2_OUTPUT_DIR}/${SAMPLE_NAME}_variants_filtered.vcf.gz

rm -f ${MUTECT2_OUTPUT_DIR}/${SAMPLE_NAME}_variants.vcf.gz  
rm -f ${MUTECT2_OUTPUT_DIR}/${SAMPLE_NAME}_variants.vcf.gz.tbi
rm -f ${MUTECT2_OUTPUT_DIR}/${SAMPLE_NAME}_variants.vcf.gz.stats

if [[ ${KEEP_ALL} -eq 0 ]]; then
    log "Selecting only variants that passed filtering..."
    run_cmd gatk SelectVariants \
        -V ${MUTECT2_OUTPUT_DIR}/${SAMPLE_NAME}_variants_filtered.vcf.gz \
        --exclude-filtered true \
        -O ${MUTECT2_OUTPUT_DIR}/${SAMPLE_NAME}_variants_filtered_passed.vcf.gz
fi

#### DOWNLOADING FUNCOTATOR DATA SOURCES #####
log "Downloading variant annotation data..."
run_cmd gatk FuncotatorDataSourceDownloader \
    --somatic \
    --validate-integrity \
    --extract-after-download \
    --${GENOME_VERS}

#### ANNOTATING DATA WITH FUNCOTATOR ####
log "Annotating variants..."

# setting variables to handle different genome versions and whether or not passing variants were selected
FUNCOTATOR_DATA_PATH="./funcotator_dataSources.v1.8.${GENOME_VERS}.20230908s"

if [[ ${KEEP_ALL} -eq 0 ]]; then 
    VCF_IN=${MUTECT2_OUTPUT_DIR}/${SAMPLE_NAME}_variants_filtered_passed.vcf.gz
    SUFFIX="passed"
else 
    VCF_IN=${MUTECT2_OUTPUT_DIR}/${SAMPLE_NAME}_variants_filtered.vcf.gz
    SUFFIX="all"
fi

OUT_VCF=${MUTECT2_OUTPUT_DIR}/${SAMPLE_NAME}_variants_filtered_${SUFFIX}_annotated.vcf.gz

run_cmd gatk Funcotator \
    --variant $VCF_IN \
    --reference $REFERENCE \
    --ref-version $GENOME_VERS \
    --data-sources-path $FUNCOTATOR_DATA_PATH \
    --output $OUT_VCF \
    --output-file-format VCF

#### CONVERTING VCF TO TABLE ####
log "Converting VCF to table..."

run_cmd gatk VariantsToTable \
    -V $OUT_VCF \
    -F CHROM -F POS -F REF -F ALT -F TYPE -F DP -F FUNCOTATION -GF AF -GF GT -GF FAD \
    -O ${MUTECT2_OUTPUT_DIR}/${SAMPLE_NAME}_variants_filtered_${SUFFIX}_annotated_table.tsv

}