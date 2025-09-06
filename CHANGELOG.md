# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2025-09-06
### Added
- Support for germline variant calling using HaplotypeCaller.
- Included hard-filtering for germline variants. 
- Updated documentation for germline analysis mode. 
  
### Changed
- Changed from using af-only-gnomad.hg38.vcf.gz resource to small_exac_common_3.hg38.vcf.gz for GetPilupSummaries. 
- Outputted total variants and passed variants VCF files.
- Changed output directory organization. 

### Removed
- Variant annotation with Funcotator.
- Output of VCF into tabular format.

## [1.0.0] - 2025-09-01
### Added
- Initial release.
- Somatic variant calling (tumor-only and tumor-normal modes).
- Automated preprocessing pipeline.
- Automatic reference database management.