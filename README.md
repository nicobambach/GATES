# GATES: GATK Automated Tool for Exome Sequencing
Version 1.1.0

GATES is a package for fully automated single-sample whole-exome sequencing (WES) analysis. It implements GATK Best Practices through a simplified command-line interface, supporting both somatic and germline variant discovery. Built for researchers without extensive computational expertise, GATES is lightweight and runs end-to-end on a standard laptop. From raw FASTQ files to filtered variants, GATES makes WES analysis accessible, reproducible, and reliable.

## Package Overview

### Sample Preprocessing
GATES supports comprehensive WES sample preprocessing with a single command and minimal user inputs, allowing users to easily go from raw FASTQ files to analysis-ready BAM files. These files can be visualized in IGV and used as direct inputs for downstream variant calling. GATES preprocessing includes: 
   - Read alignment utilizing BWA-MEM
   - Duplicate marking and BAM file sorting
   - Base quality score recalibration

### Variant Calling
GATES supports comprehensive WES variant calling and filtering with a single command, allowing users to easily go from preprocessed BAM files to variants in the VCF format. Users can call variants on preprocessed BAM files in different modes, depending on experimental question and sample availability. GATES supports:
1. **Tumor-Only Somatic Variant Calling**
   - Somatic variant calling on non-paired tumor sample using Mutect2
   - Automatically downloads and utilizes public panel of normals to identify technical artifacts
   - Automatically downloads and utilizes germline resource to exclude common germline variants
   - Estimates contamination and read orientation bias modeling for variant filtering
  
2. **Tumor-Normal Somatic Variant Calling**
   - Somatic variant calling with tumor sample and paired-normal using Mutect2
   - Automatically downloads and utilizes public panel of normals to identify technical artifacts
   - Automatically downloads and utilizes germline resource to exclude common germline variants
   - Estimates contamination and read orientation bias modeling for variant filtering

3. **Germline Variant Calling**
   - Germline variant calling using HaplotypeCaller
   - Variant filtering using GATK-recommended hard filtering thresholds

GATES automatically manages all reference databases and supporting files, requiring only input FASTQ files, a reference genome file, and capture intervals to run. The package uses conda for dependency management, ensuring reproducible analyses across different systems.

## Installation

### Prerequisites
- [Conda](https://docs.conda.io/en/latest/miniconda.html)
- Linux or macOS system
- At least 16 GB RAM recommended
- Sufficient disk space for reference files and output files (~50 GB)

### Install from GitHub

1. **Clone this repository**:
```bash
git clone https://github.com/nicobambach/GATES.git
cd GATES
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

Always create and `cd` into a new project directory before running GATES. For new projects using GATES, it is recommended to move the GATES-generated `supporting_files/` directory to the new project directory for decreased run-times (discussed in detail below).

### Tumor-Only Somatic Variant Analysis
```bash
# Preprocessing
gates preprocess \
    --sample-name SAMPLE_ID \
    --fastq1 sample_R1.fastq.gz \
    --fastq2 sample_R2.fastq.gz \
    --reference hg38.fa \
    --intervals capture_regions.bed \
    --threads 8

# Variant Calling
gates call \
    --sample-name SAMPLE_ID \
    --tumor-bam preprocessing/mapped_reads/SAMPLE_ID_recal.bam \
    --reference hg38.fa \
    --intervals capture_regions.bed \
    --mode tumor-only \
    --threads 8
```

### Tumor-Normal Somatic Variant Analysis
```bash
# Preprocessing
gates preprocess \
    --sample-name SAMPLE_ID_NORMAL \
    --fastq1 normal_R1.fastq.gz \
    --fastq2 normal_R2.fastq.gz \
    --reference hg38.fa \
    --intervals capture_regions.bed \
    --threads 8

gates preprocess \
    --sample-name SAMPLE_ID_TUMOR \
    --fastq1 tumor_R1.fastq.gz \
    --fastq2 tumor_R2.fastq.gz \
    --reference hg38.fa \
    --intervals capture_regions.bed \
    --threads 8

# Variant Calling
gates call \
    --sample-name SAMPLE_ID \
    --tumor-bam preprocessing/mapped_reads/SAMPLE_ID_TUMOR_recal.bam \
    --normal-bam preprocessing/mapped_reads/SAMPLE_ID_NORMAL_recal.bam \
    --reference hg38.fa \
    --intervals capture_regions.bed \
    --mode tumor-normal \
    --threads 8
```

### Germline Variant Analysis
```bash
# Preprocessing
gates preprocess \
    --sample-name SAMPLE_ID \
    --fastq1 sample_R1.fastq.gz \
    --fastq2 sample_R2.fastq.gz \
    --reference hg38.fa \
    --intervals capture_regions.bed \
    --threads 8

# Variant Calling
gates call \
    --sample-name SAMPLE_ID \
    --tumor-bam preprocessing/mapped_reads/SAMPLE_ID_recal.bam \
    --reference hg38.fa \
    --intervals capture_regions.bed \
    --mode germline \
    --threads 8
```
> Note: For germline analysis, use `--tumor-bam` to specify your sample BAM file. This parameter name is used across all modes for consistency.

### Command Reference

#### `gates preprocess`
Runs the preprocessing pipeline including alignment, duplicate marking, and BQSR.
```
Arguments: 
    -s, --sample-name <string>      sample identifier [REQUIRED]
        --fastq1 <path>             forward FASTQ file [REQUIRED]
        --fastq2 <path>             reverse FASTQ file [REQUIRED]
    -r, --reference <path>          reference hg38/GRCh38 FASTA file [REQUIRED]
    -i, --intervals <path>          BED/VCF/.interval_list/.list/.intervals file specifying exon capture intervals [REQUIRED]
    -t, --threads <integer>         threads to use in programs that support multithreading [OPTIONAL] [1]    

Options: 
    -v, --verbose                   display tool outputs
    -h, --help                      show help message  
```

#### `gates call`
Runs the calling pipeline including variant identification and filtering. 
```
Arguments: 
    -s, --sample-name <string>      sample identifier [REQUIRED]
        --tumor-bam <path>          preprocessed tumor BAM file [REQUIRED]
        --normal-bam <path>         preprocessed normal BAM file [REQUIRED if --mode tumor-normal]
    -r, --reference <path>          reference FASTA file [REQUIRED]
    -m, --mode <string>             mode to run pipeline. possible values: {tumor-only, tumor-normal, germline} [REQUIRED]
    -i, --intervals <path>          BED/VCF/.interval_list/.list/.intervals file specifying exon capture intervals [REQUIRED]
    -t, --threads <integer>         threads to use in programs that support multithreading [OPTIONAL] [1]    

Options: 
    -v, --verbose                   display tool outputs
    -h, --help                      show help message   
```

## Input Files

### Preprocessing Input
1. Paired-end FASTQ files (required)
2. Reference hg38/GRCh38 FASTA file (required)
3. Capture interval file: BED file or interval list defining exon capture regions (required)

### Variant Calling Input
1. Preprocessed BAM file(s) (required) (produced by `gates preprocess`)
2. Reference hg38/GRCh38 FASTA file (required)
3. Capture interval file: BED file or interval list defining exon capture regions (required)

### Automatically Downloaded Resources
GATES automatically downloads and organizes the following reference databases that are used during preprocessing and variant calling:
1. Mills_and_1000G_gold_standard.indels.hg38.vcf.gz(.tbi) (downloaded and used during `gates preprocess`)
2. Homo_sapiens_assembly38.known_indels.vcf.gz(.tbi) (downloaded and used during `gates preprocess`)
3. Homo_sapiens_assembly38.dbsnp138.vcf(.idx) (downloaded and used during `gates preprocess`)
4. 1000g_pon.hg38.vcf.gz(.tbi) (downloaded and used during `gates call`)
5. af-only-gnomad.hg38.vcf.gz(.tbi) (downloaded and used during `gates call`)
6. small_exac_common_3.hg38.vcf.gz(.tbi) (downloaded and used during `gates call`)

These files are stored in the `supporting_files/` folder in the project directory. To speed up subsequent GATES runs, keep these files in your project directory. GATES will automatically detect existing files and only re-download when necessary. You can move the entire `supporting_files/` folder to new project directories to avoid re-downloading these large files. 

## Output Files

### Preprocessing Outputs
After running `gates preprocess` the project directory is organized as follows: 
```
project_directory/
├── preprocessing/
│   ├── mapped_reads/
│   │   ├── SAMPLE_NAME_recal.bam
│   │   └── SAMPLE_NAME_recal.bai
│   └── bqsr_output/
├── supporting_files/
│   └── preprocessing_resources/
└── YYYY-MM-DD_HH-MM-SS_gates.log 
```
`gates preprocess` generates a full preprocessed BAM file (aligned, sorted, duplicate-marked, and BQSR-adjusted) that is now ready for variant calling. The outputted BAM file(s) from preprocessing runs are used as inputs for variant calling using `gates call`. It is important to keep the `*.bai` file in the same directory as the BAM file as it is needed for many downstream tools. These index files are also necessary for viewing the preprocessed BAM in IGV. 

### Variant Calling Outputs
After running `gates call` the project directory is organized as follows:
```
project_directory/
├── preprocessing/
├── MODE/
│   └── variants/
│       ├── SAMPLE_NAME_all_germline_variants.vcf.gz
│       └── SAMPLE_NAME_passed_germline_variants.vcf.gz
├── supporting_files/
│   ├── preprocessing_resources/
│   └── calling_resources/
└── YYYY-MM-DD_HH-MM-SS_gates.log 
```
`gates call` generates files containing variants called by either Mutect2 (somatic) or HaplotypeCaller (germline) in the VCF file format. The output folder is named based on analysis mode (`germline/`, `tumor_only_somatic/`, or `tumor_normal_somatic/`) and VCF file names include whether called variants are germline or somatic. 

Variants are automatically filtered based on GATK-recommended parameters for germline and somatic variant discovery. `SAMPLE_NAME_passed_*_variants.vcf.gz` contains only variants that passed filtering, representing high-confidence calls. `SAMPLE_NAME_all_*_variants.vcf.gz` includes both passing and non-passing variants with the `FILTER` field annotated based on which filtering criteria each variant did not pass. This file may be useful for quality control, however, it is recommended to use only passing variants for downstream analysis. 

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
tail -50 YYYY-MM-DD_HH-MM-SS_gates.log
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Support

- **Issues**: Report bugs and request features on [GitHub Issues](https://github.com/nicobambach/gates/issues)
- **Documentation**: Check this README and command help (`gates --help`, `gates preprocess --help`, `gates call --help`)

## Citation

If you use GATES in your research, please cite:

```
Bambach, NE (2025). GATES: GATK Automated Tool for Exome Sequencing v1.1.0. 
GitHub. https://github.com/nicobambach/gates
```