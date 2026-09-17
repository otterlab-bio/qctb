use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};
use qctb::qc_summary;
use std::path::{Path, PathBuf};

#[derive(Clone, Copy, Debug, Eq, PartialEq, ValueEnum)]
#[value(rename_all = "lower")]
enum OutputFormat {
    Xlsx,
    Tsv,
}

impl OutputFormat {
    fn as_str(self) -> &'static str {
        match self {
            Self::Xlsx => "xlsx",
            Self::Tsv => "tsv",
        }
    }
}

#[derive(Parser)]
#[command(name = "qctb")]
#[command(about = "QC tools for bioinformatics", long_about = None)]
#[command(version)]
struct Cli {
    /// YAML configuration file. Mutually exclusive with --config-dir.
    #[arg(
        long,
        required_unless_present = "config_dir",
        conflicts_with = "config_dir"
    )]
    config: Option<String>,

    /// Otter run directory containing immutable run.yaml.
    #[arg(long, required_unless_present = "config", conflicts_with = "config")]
    config_dir: Option<String>,

    /// Output file path
    #[arg(long)]
    output: String,

    /// Output format: lowercase xlsx or tsv (default: xlsx)
    #[arg(long, value_enum, default_value = "xlsx", ignore_case = false)]
    format: OutputFormat,

    /// RNA-seq mode (use RNA-seq specific metrics and parsers)
    #[arg(long)]
    rnaseq: bool,
}

fn resolve_config_path(cli: &Cli) -> Result<PathBuf> {
    if let Some(config_path) = &cli.config {
        return Ok(PathBuf::from(config_path));
    }
    let run_directory = Path::new(
        cli.config_dir
            .as_deref()
            .expect("clap requires either --config or --config-dir"),
    );
    if !run_directory.is_dir() {
        anyhow::bail!(
            "QCTB --config-dir must name an Otter run directory: {}",
            run_directory.display()
        );
    }
    Ok(run_directory.join("run.yaml"))
}

fn workflow_uses_rnaseq_parser(mode: &str) -> bool {
    matches!(
        mode.trim().to_ascii_uppercase().as_str(),
        "RNASEQ" | "RNA-SEQ" | "RNA-PDX" | "BEAVERRNASEQPDX"
    )
}

fn resolve_report_mode(
    workflow_mode: Option<&str>,
    rnaseq: bool,
) -> Result<qc_summary::ReportMode> {
    let declared_mode = workflow_mode.map(str::trim).filter(|mode| !mode.is_empty());
    let declared_rnaseq = declared_mode.is_some_and(workflow_uses_rnaseq_parser);

    if declared_rnaseq && !rnaseq {
        anyhow::bail!(
            "Configuration declares an RNA-seq workflow; invoke qctb with --rnaseq to select RNA parsers"
        );
    }
    if rnaseq && !declared_rnaseq {
        if let Some(mode) = declared_mode {
            anyhow::bail!("--rnaseq conflicts with configured workflow mode '{mode}'");
        }
    }

    Ok(qc_summary::ReportMode::from_workflow_mode(
        declared_mode,
        rnaseq,
    ))
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let config_path = resolve_config_path(&cli)?;
    let qc_config = qc_summary::load_config(&config_path)
        .with_context(|| format!("Failed to load config from: {}", config_path.display()))?;
    let report_mode = resolve_report_mode(qc_config.workflow_mode.as_deref(), cli.rnaseq)?;

    println!(
        "Processing {} samples in {} mode...",
        qc_config.SIDs.len(),
        report_mode.as_str()
    );
    println!("Output format: {}", cli.format.as_str());
    println!("Output schema: {}", qc_summary::REPORT_SCHEMA_ID);

    if report_mode == qc_summary::ReportMode::RnaSeq {
        let summaries = qc_summary::process_all_samples_rnaseq(&qc_config)
            .with_context(|| "Failed to process samples in RNA-seq mode")?;
        println!("Successfully processed {} samples", summaries.len());
        match cli.format {
            OutputFormat::Xlsx => qc_summary::write_excel_rnaseq(&summaries, &cli.output)
                .with_context(|| format!("Failed to write Excel output to: {}", cli.output))?,
            OutputFormat::Tsv => {
                qc_summary::write_tsv_rnaseq(&summaries, Path::new(&cli.output))
                    .with_context(|| format!("Failed to write TSV output to: {}", cli.output))?
            }
        }
    } else {
        let summaries = qc_summary::process_all_samples(&qc_config)
            .with_context(|| "Failed to process samples")?;
        println!("Successfully processed {} samples", summaries.len());
        match cli.format {
            OutputFormat::Xlsx => {
                qc_summary::write_excel_standard_mode(&summaries, &cli.output, report_mode)
                    .with_context(|| format!("Failed to write Excel output to: {}", cli.output))?
            }
            OutputFormat::Tsv => {
                qc_summary::write_tsv_standard(&summaries, Path::new(&cli.output), report_mode)
                    .with_context(|| format!("Failed to write TSV output to: {}", cli.output))?
            }
        }
    }

    println!("Output written to: {}", cli.output);
    Ok(())
}

#[cfg(test)]
mod cli_tests {
    use super::*;

    fn base_arguments() -> [&'static str; 5] {
        [
            "qctb",
            "--config",
            "config.yaml",
            "--output",
            "summary.xlsx",
        ]
    }

    #[test]
    fn config_file_and_run_directory_inputs_are_mutually_exclusive() {
        let directory_cli = Cli::try_parse_from([
            "qctb",
            "--config-dir",
            "/analysis/runs/run-example",
            "--output",
            "summary.xlsx",
        ])
        .expect("canonical run directory should parse");
        assert_eq!(
            directory_cli.config_dir.as_deref(),
            Some("/analysis/runs/run-example")
        );
        assert!(directory_cli.config.is_none());

        let conflicting_cli = Cli::try_parse_from([
            "qctb",
            "--config",
            "run.yaml",
            "--config-dir",
            "/analysis/runs/run-example",
            "--output",
            "summary.xlsx",
        ]);
        assert!(conflicting_cli.is_err());
    }

    #[test]
    fn rnaseq_flag_must_match_the_declared_workflow_mode() {
        assert_eq!(
            resolve_report_mode(Some("RNASEQ"), true).unwrap(),
            qc_summary::ReportMode::RnaSeq
        );
        assert_eq!(
            resolve_report_mode(Some("BEAVERRNASEQPDX"), true).unwrap(),
            qc_summary::ReportMode::RnaSeq
        );
        assert_eq!(
            resolve_report_mode(Some("RNA-PDX"), true).unwrap(),
            qc_summary::ReportMode::RnaSeq
        );
        assert!(resolve_report_mode(Some("RNASEQ"), false).is_err());
        assert!(resolve_report_mode(Some("RRBS"), true).is_err());
        assert!(resolve_report_mode(Some("BEAVERPDX"), true).is_err());

        assert_eq!(
            resolve_report_mode(None, true).unwrap(),
            qc_summary::ReportMode::RnaSeq
        );
        assert_eq!(
            resolve_report_mode(None, false).unwrap(),
            qc_summary::ReportMode::Standard
        );
    }

    #[test]
    fn output_format_accepts_only_lowercase_xlsx_or_tsv() {
        let default_cli =
            Cli::try_parse_from(base_arguments()).expect("default format should parse");
        assert_eq!(default_cli.format, OutputFormat::Xlsx);

        let tsv_cli = Cli::try_parse_from([
            "qctb",
            "--config",
            "config.yaml",
            "--output",
            "summary.tsv",
            "--format",
            "tsv",
        ])
        .expect("lowercase tsv should parse");
        assert_eq!(tsv_cli.format, OutputFormat::Tsv);

        for invalid_format in ["json", "XLSX", "Tsv"] {
            let result = Cli::try_parse_from([
                "qctb",
                "--config",
                "config.yaml",
                "--output",
                "summary.out",
                "--format",
                invalid_format,
            ]);
            assert!(
                result.is_err(),
                "format '{invalid_format}' should be rejected"
            );
        }
    }
}
