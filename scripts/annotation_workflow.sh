#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# getting root directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# sourcing helper file to load run_cmd and log functions
source "$SCRIPT_DIR/scripts/helper.sh"

annotation_workflow() {

# setting directories
VCF_DIR="$(dirname "$(realpath "$VCF")")" 
ANNOTATION_VARIANTS_DIR="$VCF_DIR/annotated_variants"

# making necessary directories
mkdir -p "$ANNOTATION_VARIANTS_DIR"

# running VEP
run_cmd vep \
    --input_file $VCF \
    --format vcf \
    --output_file ${ANNOTATION_VARIANTS_DIR}/${SAMPLE_NAME}_all_variants_annotated.vcf \
    --vcf \
    --stats_file ${ANNOTATION_VARIANTS_DIR}/${SAMPLE_NAME}_vep_stats.html \
    --stats_html \
    --cache --dir_cache $CACHE \
    --refseq \
    --fasta $REFERENCE \
    --variant_class \
    --sift b \
    --polyphen b \
    --mane \
    --total_length \
    --biotype \
    --numbers \
    --domains \
    --symbol \
    --pick \
    --variant_class \
    --hgvs \
    --max_af \
    --af_gnomade \
    --af_1kg \
    --force_overwrite \
    --offline

# filtering VCF to only protein coding variants and by population AF threshold
run_cmd filter_vep \
    --input_file ${ANNOTATION_VARIANTS_DIR}/${SAMPLE_NAME}_all_variants_annotated.vcf \
    --format vcf \
    --output_file ${ANNOTATION_VARIANTS_DIR}/${SAMPLE_NAME}_rare_nonsyn_variants.vcf \
    --filter "(MAX_AF < ${POP_AF} or not MAX_AF) and Consequence != synonymous_variant and Consequence != intron_variant and Consequence != downstream_gene_variant and Consequence != upstream_gene_variant" \
    --force_overwrite

# converting filtered VCF to TSV
if [[ "$MODE" == "tumor-normal" ]]; then

run_cmd bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%VARIANT_CLASS\t%Consequence\t%Feature\t%EXON\t%HGVSc\t%HGVSp\t%cDNA_position\t%Protein_position\t%Amino_acids\t%STRAND\t%SIFT\t%PolyPhen[\t%SAMPLE=%AF\t%SAMPLE=%AD]\n' -HH ${ANNOTATION_VARIANTS_DIR}/${SAMPLE_NAME}_rare_nonsyn_variants.vcf -o ${ANNOTATION_VARIANTS_DIR}/${SAMPLE_NAME}_rare_nonsyn_variants.tsv

elif [[ "$MODE" == "germline" ]]; then

run_cmd bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%VARIANT_CLASS\t%Consequence\t%Feature\t%EXON\t%HGVSc\t%HGVSp\t%cDNA_position\t%Protein_position\t%Amino_acids\t%STRAND\t%SIFT\t%PolyPhen\t%INFO/AF[\t%AD]\n' -HH ${ANNOTATION_VARIANTS_DIR}/${SAMPLE_NAME}_rare_nonsyn_variants.vcf -o ${ANNOTATION_VARIANTS_DIR}/${SAMPLE_NAME}_rare_nonsyn_variants.tsv


elif [[ "$MODE" == "tumor-only" ]]; then

run_cmd bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%SYMBOL\t%VARIANT_CLASS\t%Consequence\t%Feature\t%EXON\t%HGVSc\t%HGVSp\t%cDNA_position\t%Protein_position\t%Amino_acids\t%STRAND\t%SIFT\t%PolyPhen[\t%AF\t%AD]\n' -HH ${ANNOTATION_VARIANTS_DIR}/${SAMPLE_NAME}_rare_nonsyn_variants.vcf -o ${ANNOTATION_VARIANTS_DIR}/${SAMPLE_NAME}_rare_nonsyn_variants.tsv

fi

PASSED_VARIANTS=$(bcftools view -H ${ANNOTATION_VARIANTS_DIR}/${SAMPLE_NAME}_rare_nonsyn_variants.vcf| wc -l)

log "Finished annotating variants for $SAMPLE_NAME"
log "Number of rare, non-synonumous variants: $PASSED_VARIANTS"
log "Variant annotation output files:"
log "  - ${SAMPLE_NAME}_all_variants_annotated.vcf (VCF of all variants annotated)"
log "  - ${SAMPLE_NAME}_rare_nonsyn_variants.vcf (annotated VCF with rare, non-synonymous variants)"
log "  - ${SAMPLE_NAME}_rare_nonsyn_variants.tsv (TSV of rare, non-synonymous variants)"

}








