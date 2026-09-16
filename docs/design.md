# qctb Architecture Design

## Data Flow

```text
config.yaml
  -> dual-format config loader + SID validation
  -> fastqcx / Trim / Bismark / Qualimap / STAR / methx parsers
  -> QCSummary or QCSummaryRNA
  -> typed ReportRow using qctb.report/1.0.0
  -> atomic Excel or TSV publication
```

## Native Report Contract

`src/qc_summary/schema.rs` is the single source of truth for both writers. `ReportColumnSpec` freezes each column's machine key, Excel header, type, decimal places, order, and `N/A` missing representation. `ReportCell` enforces that values match the declared type before publication.

RRBS, WGBS, and PDX v1 share the 42-column BS report layout. Their mode identity remains distinct in metadata and golden tests. RNA-seq uses a separate 22-column layout. Adding PDX host-specific fields or otherwise changing a layout requires a schema version bump.

TSV embeds:

```text
# qctb_schema=qctb.report/1.0.0
# qctb_mode=<mode>
```

Excel keeps `Report` first for normal spreadsheet use and adds `qctb_metadata` with schema and column definitions.

## Parser Policy

Required fields must occur exactly once. Counts use integer parsing, percentages and other decimals must be finite and in their declared range, and derived metrics must agree with source counts within display-rounding tolerance. Contradictory reports fail closed instead of producing partial summaries.

The `fastqcx` parser consumes the FastQC-compatible `fastqc_data.txt` format and the embedded `Seqkit Statistics` module. Its module, public symbols, output-directory suffixes, and test fixtures use only the `fastqcx` product name.

The `methx` integration consumes Methrix scientific workbooks using only the declared sheets and headers:

- `Sheet1`: coverage columns from `methx`.
- `ChIPseeker_By_Sample`: `sample`, `covered_cpgs`, and the contracted Promoter/Exon/Intron/Intergenic count and percent pairs.

Extra annotation columns may be retained by the parser, but the v1 report emits only contracted metrics. Missing contracted metrics fail. Existing Methrix report labels, schema keys, workbook names, and the `methrixh5` data directory remain unchanged.

## Publication

`write_atomically` creates the temporary output in the destination directory, writes the complete workbook or TSV, synchronizes it, and atomically replaces the target. A write failure leaves an existing target untouched. Parent directories are created before staging.

## Testing

- Parser unit tests cover missing, duplicate, invalid-range, overflow, and inconsistent fields.
- Methrix workbooks are generated at test runtime, so a clean clone does not depend on ignored binary fixtures.
- RRBS, WGBS, RNA-seq, and PDX each have a byte-stable TSV golden.
- Excel integration tests read back the metadata and report sheets.
