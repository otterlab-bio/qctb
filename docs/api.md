# qctb CLI API Reference

## Overview

qctb provides bioinformatics QC tools implemented in Rust for high performance.

## Command

Generate QC summary report from multiple tool outputs.

```bash
qctb [OPTIONS]
```

**Options**:
- `--config <FILE>` - YAML configuration file, including canonical immutable `otter.run/v1` snapshots (required unless `--config-dir` is supplied)
- `--config-dir <DIRECTORY>` - Otter run directory containing `run.yaml` (required unless `--config` is supplied)
- `--output <FILE>` - Output file path (required)
- `--format <fmt>` - Output format: xlsx or tsv (default: xlsx)
- `--rnaseq` - Select RNA-seq parsers; required when the configuration declares RNA-seq or RNA-PDX, and rejected for other explicit modes

**Standard Mode**: Parses FastQC-compatible `fastqcx` Seqkit Statistics, Trim Galore, Bismark, Qualimap, and optional Methrix QC outputs produced by `methx`.
**RNA-seq Mode**: Parses FastQC-compatible `fastqcx` Seqkit Statistics, Trim Galore, and STAR outputs.

**Exit Codes**:
- `0` - Success
- `1` - Error (check stderr for details)

**Examples**:
```bash
# Standard QC summary (BS-seq)
qctb \
  --config config/config.yaml \
  --output qc_summary.xlsx

# Standard QC summary (TSV format)
qctb \
  --config config/config.yaml \
  --output qc_summary.tsv \
  --format tsv

# RNA-seq QC summary
qctb \
  --config config/config.yaml \
  --output qc_summary_rnaseq.xlsx \
  --rnaseq
```

## Global Options

- `-h, --help` - Print help information
- `-V, --version` - Print version information

## Configuration File Format

The YAML config file must contain:

```yaml
SIDs:
  - sample1
  - sample2
qcDir: "/path/to/qc/output"
trimDir: "/path/to/trim/output"
bsmapDir: "/path/to/bsmap"
outDir_mCall: "/path/to/mcall"
# Optional fields
graft: "human"  # or other species
qcdir_before: "/path/to/before/qc"
qcdir_after: "/path/to/after/qc"
```

## Output Formats

All outputs implement `qctb.report/1.0.0`. Schema changes require an explicit version bump.

### Excel Output
- `Report` is the first worksheet and contains the typed report table.
- `qctb_metadata` records schema name/version, mode, missing-value marker, and every column's key, header, type, and decimal places.
- RRBS, WGBS, and PDX v1 use the 42-column BS schema; RNA-seq uses 22 columns.

### TSV Output
- UTF-8 tab-separated output.
- Line 1: `# qctb_schema=qctb.report/1.0.0`.
- Line 2: `# qctb_mode=<RRBS|WGBS|RNA-seq|PDX|standard>`.
- Line 3 is the stable machine column header.
- Missing optional values are exactly `N/A`.
- Consumers should treat lines beginning with `#` as metadata comments.
