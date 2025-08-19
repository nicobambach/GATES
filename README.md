# GATES: GATK Automated Tool for Exome Sequencing

GATES is a lightweight, fully automated pipeline for single-sample somatic whole exome sequencing (WES) analysis. It implements both the GATK data pre-processing and somatic short variant discovery best practices workflow with two simple commands. Built for wet-lab scientists with limited computational resources, GATES runs end-to-end on a standard laptop so basic WES analysis is simple, reproducible, and accessible to everyone.

## Pipeline Overview

### Preprocessing
1. **Read Alignment**: BWA-MEM alignment to hg38/GRCh38 reference
2. **Duplicate Marking**: GATK MarkDuplicatesSpark
3. **Base Quality Score Recalibration (BQSR)**: Recalibration of base quality scores based on known variant sites
4. **Quality Assessment**: Comprehensive plots comparing base qualities before and after BQSR

### Variant Calling Pipeline
1. **Somatic Variant Calling**: GATK Mutect2 with publicly available Panel of Normals
2. **Contamination Estimation**: Cross-sample contamination analysis
3. **Orientation Bias Modeling**: Read orientation artifact detection
4. **Variant Filtering**: GATK FilterMutectCalls with multiple filters
5. **Variant Annotation**: Variant annotation using GATK Funcotator 
6. **Output Generation**: Both VCF and tab-delimited table formats

## Installation

### Prerequisites
- [Conda](https://docs.conda.io/en/latest/miniconda.html)
- Linux or macOS system
- At least 16GB RAM recommended
- Sufficient disk space for reference files and output files (~50GB)

### Install from GitHub

1. **Clone this repository**:
```bash
git clone https://github.com/nicobambach/gates.git
cd gates
```

2. **Create and activate conda environment**:
```bash
conda env create -f environment.yaml
conda activate gates
```

3. **Make scripts executable and add to PATH**:
```bash
chmod +x bin/gates scripts/*.sh
export PATH="$PWD/bin:$PATH"
```

4. **Test installation**:
```bash
gates --version
gates --help
```

## Usage

### Tumor-Only Analysis
```bash
# Preprocessing
gates preprocess \
    --sample-name SAMPLE_ID \
    --tumor-fq1 tumor_R1.fastq.gz \
    --tumor-fq2 tumor_R2.fastq.gz \
    --reference hg38.fa \
    --intervals capture_regions.bed \
    --mode tumor-only \
    --threads 8

# Variant Calling
gates call \
    --sample-name SAMPLE_ID \
    --tumor-bam mapped_reads/SAMPLE_ID_sorted_dedup_recal.bam \
    --reference hg38.fa \
    --intervals capture_regions.bed \
    --mode tumor-only
    --threads 8
```

### Tumor-Normal Analysis
```bash
# Preprocessing
gates preprocess \
    --sample-name SAMPLE_ID \
    --tumor-fq1 tumor_R1.fastq.gz \
    --tumor-fq2 tumor_R2.fastq.gz \
    --normal-fq1 normal_R1.fastq.gz \
    --normal-fq2 normal_R2.fastq.gz \
    --reference hg38.fa \
    --intervals capture_regions.bed \
    --mode tumor-normal \
    --threads 8

# Variant Calling
gates call \
    --sample-name SAMPLE_ID \
    --tumor-bam mapped_reads/SAMPLE_ID_tumor_sorted_dedup_recal.bam \
    --normal-bam mapped_reads/SAMPLE_ID_normal_sorted_dedup_recal.bam \
    --reference hg38.fa \
    --intervals capture_regions.bed \
    --mode tumor-normal \
    --threads 8
```

### Command Reference

#### `gates preprocess`
Runs the preprocessing pipeline including alignment, duplicate marking, and BQSR.

**Arguments:**
- `--sample-name`: Sample identifier
- `--tumor-fq1/2`: Paired-end tumor FASTQ files
- `--reference`: hg38/GRCh38 reference FASTA file
- `--intervals`: BED file or interval list for capture regions
- `--mode`: `tumor-only` or `tumor-normal`

**Tumor-Normal Mode Additional Requirements:**
- `--normal-fq1/2`: Paired-end normal FASTQ files

**Optional Arguments:**
- `--threads`: Number of threads (default: 1)
- `--verbose`: Show tool outputs in terminal window (note: all tool outputs are captured in time-stamped log files)
- `--help`: Show usage information

**Arguments:**
    `-s, --sample-name` <string>      sample identifier [REQUIRED]
        `--tumor-fq1` <path>          tumor forward FASTQ file [REQUIRED]
        `--tumor-fq2` <path>          tumor reverse FASTQ file [REQUIRED]
        `--normal-fq1` <path>         normal forward FASTQ file [REQUIRED if --mode tumor-normal]
        `--normal-fq2` <path>         normal reverse FASTQ file [REQUIRED if --mode tumor-normal]
    `-r, --reference` <path>          reference hg38/GRCh38 FASTA file [REQUIRED]
    `-m, --mode` <string>             mode to run preprocessing. possible values: {tumor-only, tumor-normal} [REQUIRED]
    `-i, --intervals` <path>          BED/VCF/.interval_list/.list/.intervals file specifying exon capture intervals [REQUIRED]
    `-t, --threads` <integer>         threads to used in programs that support multithreading [OPTIONAL] [1]    

Options: 
`-v, --verbose`                     display tool outputs
`-h, --help`                         show help message                         

#### `gates call`
Runs the variant calling pipeline including Mutect2, filtering, and annotation.

**Required Arguments:**
- `--sample-name`: Sample identifier
- `--tumor-bam`: Preprocessed tumor BAM file (outputted from `gates preprocess`)
- `--reference`: Reference FASTA file
- `--intervals`: BED file or interval list for capture regions
- `--mode`: `tumor-only` or `tumor-normal`

**Tumor-Normal Mode Additional Requirements:**
- `--normal-bam`: Preprocessed normal BAM file (outputted from `gates preprocess`)

**Optional Arguments:**
- `--keep-all`: Retain all variants, including those that don't pass filtering (not recommended)
- `--threads`: Number of threads (default: 1)
- `--verbose`: Show detailed tool outputs in terminal window (note: all tool outputs are captured in time-stamped log files)
- `--help`: Show usage information

## Input Files

### Preprocessing Input
1. Paired-end tumor FASTQ files (required)
2. Paired-end normal FASTQ files (required for tumor-normal mode)
3. Reference hg38/GRCh38 FASTA file (required)
4. Capture interval file: BED file or interval list defining exon capture regions (required)

### Variant Calling Input
1. Processed tumor BAM file (required) (outputted from `gates preprocess`)
2. Procesed normal BAM file (required for tumor-normal mode) (outputted from `gates preprocess`)
3. Reference hg38/GRCh38 FASTA file (required)
4. Capture interval file: BED file or interval list defining exon capture regions (required)

### Automatically Downloaded Resources
GATES automatically downloads and organizes the following reference databases that are used during preprocessing and variant calling:
1. Mills_and_1000G_gold_standard.indels.hg38.vcf.gz(.tbi) (downloaded and used during `gates preprocess`)
2. Homo_sapiens_assembly38.known_indels.vcf.gz(.tbi) (downloaded and used during `gates preprocess`)
3. Homo_sapiens_assembly38.dbsnp138.vcf(.idx) (downloaded and used during `gates preprocess`)
4. 1000g_pon.hg38.vcf.gz(.tbi) (downloaded and used during `gates call`)
5. af-only-gnomad.hg38.vcf.gz(.tbi) (downloaded and used during `gates call`)

These files are housed in `project_directory/bqsr_training` and `project_directory/mutect2_supporting_files`. To speed up subsequent GATES runs, do not delete these files and keep them in your project directory. GATES will automatically detect if these directories and files exist and only redownload when necessary.

## Output Files

### Preprocessing Outputs
After running `gates preprocess` the project directory is organized as follows: 
```
project_directory/
├── mapped_reads/
│   ├── *_sorted_dedup_recal.bam
│   └── *.bai
├── bqsr_output/
│   ├── *_recalibration_table.table
│   └── *_recalibration_plots.pdf
├── bqsr_training/
│   ├── Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
│   ├── Homo_sapiens_assembly38.known_indels.vcf.gz                   
│   ├── Homo_sapiens_assembly38.dbsnp138.vcf
│   ├── *.tbi
│   └── *.idx
└── *.log
```
The most important output files are `*_sorted_dedup_recal.bam`. These fully preprocessed BAMs (sorted, duplicate-marked, and BQSR-adjusted) are now ready for variant calling. Depending on whether `gates` is run in `tumor-only` or `tumor-normal` mode, there will be one or two BAM files. These files will be used as inputs for `gates call`. It is important to keep the `*.bai` files in this same directory as the BAM files as these are needed for downstream tools to properly read the BAM files. These index files are also necessary for viewing the aligned BAM files in IGV. 

### Variant Calling Outputs
After running `gates call` the project directory is organized as follows:
```
project_directory/
├── mutect2_output/
│   ├── *_variants_filtered_passed_annotated.vcf.gz
│   ├── *_variants_filtered_passed_annotated_table.tsv
│   ├── *_filtering_stats.txt
│   └── *.tbi
├── mutect2_filtering_data/
│   ├── *_contamination_table.table
│   ├── *_tumor_segmentation_table.table
│   └── *_read_orientation_model.tar.gz
├── mutect2_supporting_files/
│   ├── 1000g_pon.hg38.vcf.gz
│   ├── af-only-gnomad.hg38.vcf.gz
│   ├── *.tbi
│   └── *.idx
├── funcotator_dataSources.v1.8.hg38.*/
│   └── (annotation databases)
└── *.log
```
The most important output files are `*_variants_filtered_passed_annotated.vcf.gz` and `*_variants_filtered_passed_annotated_table.tsv`. These contain the final high-confidence somatic variants in both VCF format (for downstream analysis) and tab-delimited format (for easy viewing in Excel). The VCF file includes all variant annotations from Funcotator, while the TSV table provides the same information in a more accessible format. The `*_filtering_stats.txt` file provides a summary of how many variants passed or failed each filtering step, which is useful for quality assessment. If the `--keep-all` flag is used, all variants, even those that did not pass filtering, are kept and annotated. In this case, the `*_variants_filtered_passed_annotated.vcf.gz` and `*_variants_filtered_passed_annotated_table.tsv` are instead named `*_variants_filtered_all_annotated.vcf.gz` and `*_variants_filtered_all_annotated_table.tsv`, respectively. Using the `--keep-all` flag is not recommended and will cause increased run times. 

## Dependencies

All dependencies are automatically managed through conda via the environment.yaml file. It is essential to create the gates conda environment using this .yaml file. 

## Troubleshooting

### Common Issues

1. **Memory errors**: Increase available RAM or reduce thread count
2. **Download failures**: Check internet connection and retry
3. **Permission errors**: Ensure scripts are executable (`chmod +x`)
4. **Path issues**: Verify conda environment activation

### Log Files
GATES creates timestamped log files for each run. Check these for detailed error messages:
```bash
ls *.log
tail -50 [latest_log_file]
```

## Citation

If you use GATES in your research, please cite:

```
Bambach, NE (2025). GATES: GATK Automated Tool for Exome Sequencing. 
GitHub. https://github.com/nicobambach/gates
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Support

- **Issues**: Report bugs and request features on [GitHub Issues](https://github.com/nicobambach/gates/issues)
- **Documentation**: Check this README and command help (`gates --help`, `gates preprocess --help`, `gates call --help`)

---