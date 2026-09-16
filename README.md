# qctb

**A Rust QC summary reporter for Otter sequencing workflows.**

`qctb` parses FastQC-compatible `fastqc_data.txt` files from `fastqcx`, reads QC inputs produced by `methx`, and writes versioned Excel or TSV summaries for bisulfite, RNA-seq, and PDX workflows.

## Proof

Both modes write a report that identifies its own schema and mode before any row is parsed. Real
output from the tracked fixtures, truncated at the ratio column:

```text
$ qctb --config config.yaml --output qc_summary.tsv --format tsv
Output schema: qctb.report/1.0.0
Successfully processed 1 samples
Output written to: qc_summary.tsv

$ head -4 qc_summary.tsv
# qctb_schema=qctb.report/1.0.0
# qctb_mode=PDX
sample_id	reads_raw	bases_raw	reads_clean	bases_clean	clean_data_ratio	...
F6703_372760	48130550	7219582500	47108944	6614755204	0.9162	...
```

The two leading comment lines are the contract: a consumer can assert the schema version and the
mode instead of inferring them from the columns. RNA-seq mode writes `# qctb_mode=RNA-seq` in the
same position.

## Quick start

Build from source:

```bash
git clone https://github.com/otterlab-bio/qctb.git
cd qctb
cargo build --release
```

Generate a standard report:

```bash
qctb --config config.yaml --output qc_summary.xlsx
```

Generate an RNA-seq report:

```bash
qctb --config config.yaml --output qc_summary.xlsx --rnaseq
```

Use an immutable Otter run snapshot or its run directory:

```bash
qctb \
  --config /analysis/runs/<run-id>/run.yaml \
  --output qc_summary.xlsx

qctb \
  --config-dir /analysis/runs/<run-id> \
  --output qc_summary.tsv \
  --format tsv
```

`--config-dir` resolves only `<run-directory>/run.yaml` and requires the canonical immutable run contract. Legacy/nested YAML remains supported for compatibility callers through `--config`.

## CLI contract

| Option | Meaning |
| --- | --- |
| `--config` | Configuration or immutable run snapshot; mutually exclusive with `--config-dir`. |
| `--config-dir` | Otter run directory containing `run.yaml`. |
| `--output` | Required Excel or TSV output path. |
| `--format` | Lowercase `xlsx` or `tsv`; default is `xlsx`. |
| `--rnaseq` | Select RNA-seq-specific parsers and metrics. |

## Output

The report schema is emitted in the command output and encoded by the report writer. Excel is the default for workflow delivery; TSV is available for inspection and downstream automation.

QCTB is a reporting layer, not a scientific parity comparator. Compare domain metrics with the owning workflow evidence and preserve the input/reference identity.

## Development

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features --locked -- -D warnings
cargo test --all-targets --all-features --locked
cargo build --release
```

## Continuous integration

`.github/workflows/ci.yml` runs on `push`, `pull_request`, and manual dispatch. It runs
`cargo fmt --check`, Clippy with warnings denied, the full test suite, and a release
build. It then produces real report artifacts from the tracked fixtures: a standard PDX
TSV and XLSX plus an RNA-seq TSV, and asserts the `qctb.report/1.0.0` schema header and
mode for each. All logs and reports are uploaded as the `qctb-report-evidence` artifact,
including on failure.

## License and repository

MIT · [otterlab-bio/qctb](https://github.com/otterlab-bio/qctb)
