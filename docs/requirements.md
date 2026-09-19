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
- No cross-implementation speed ratio is part of the release contract.
- Any performance claim must name the fixture, version, hardware, command, and retained benchmark evidence.

### Compatibility
- Output uses the versioned native schema `qctb.report/1.0.0`
  (`src/qc_summary/schema.rs`); it is not a byte-for-byte reproduction of the
  legacy R output. Parity with the legacy R report is asserted only for the
  columns and fixtures covered by the golden tests (`tests/golden/`).
- Must accept the documented input report formats (FastQC-compatible
  `fastqc_data.txt`, Trim Galore, Bismark, Qualimap, STAR, optional methx).
- Must accept legacy `--config` YAML and canonical `--config-dir` run
  directories (`run.yaml` with `run.immutable: true` and absolute
  `paths.work`).

### Reliability
- Optional report files may be absent; a present report with missing,
  duplicate, or out-of-range required fields fails closed rather than being
  skipped.
- Must provide clear error messages
- Must validate input before processing

## Technical Constraints

### Dependencies
- clap 4.5 (CLI)
- serde/serde_yaml (YAML parsing)
- rust_xlsxwriter (Excel output)
- calamine (Excel input)
- tempfile (atomic output staging)
- anyhow (error handling)

### Platform
- Primary: Linux x86_64
- Build: stable Rust with the committed `Cargo.lock`

## Integration Requirements

### CLI Interface
```bash
qctb --config <config.yaml> --output <output.xlsx>
# The configuration must declare RNA-seq or RNA-PDX.
qctb --config <config.yaml> --output <output.xlsx> --rnaseq
qctb --config <config.yaml> --output <output.tsv> --format tsv
```

### Snakemake Integration
- Must be callable from Snakemake rules
- Must return appropriate exit codes
- Must log progress to stdout/stderr
