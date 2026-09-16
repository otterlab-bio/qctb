# qctb Technical Requirements

## Overview
Rust implementation of QC summary reporting for performance and maintainability.

## Functional Requirements

### QC Summary (Standard Mode)
- Parse YAML configuration file
- Read and parse 4 types of QC outputs:
  - FastQC-compatible `fastqc_data.txt` files produced by `fastqcx`
  - trim galore reports
  - bismark alignment results
  - qualimap reports
- Aggregate metrics across all samples
- Generate Excel/TSV summary report

### QC Summary (RNA-seq Mode)
- Parse same YAML configuration
- Read RNA-seq specific outputs:
  - FastQC-compatible `fastqc_data.txt` files produced by `fastqcx`
  - trim galore reports
  - STAR alignment results
- Aggregate RNA-seq metrics
- Generate Excel/TSV summary report

## Non-Functional Requirements

### Performance
- QC Summary: 5-6x faster than R version (30-60s → 5-10s for 10 samples)

### Compatibility
- Output must match R version exactly (same Excel format, same numbers)
- Must accept same input file formats
- Must work with existing config.yaml structure

### Reliability
- Must handle missing files gracefully
- Must provide clear error messages
- Must validate input before processing

## Technical Constraints

### Dependencies
- clap 4.5 (CLI)
- serde/serde_yaml (YAML parsing)
- rust_xlsxwriter (Excel output)
- rayon (Parallel processing)
- anyhow/thiserror (Error handling)

### Platform
- Primary: Linux x86_64
- Build: Rust 1.70+ with stable toolchain

## Integration Requirements

### CLI Interface
```bash
qctb --config <config.yaml> --output <output.xlsx>
qctb --config <config.yaml> --output <output.xlsx> --rnaseq
qctb --config <config.yaml> --output <output.tsv> --format tsv
```

### Snakemake Integration
- Must be callable from Snakemake rules
- Must return appropriate exit codes
- Must log progress to stdout/stderr
