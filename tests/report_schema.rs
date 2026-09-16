use anyhow::Result;
use calamine::{open_workbook_auto, Reader};
use qctb::qc_summary::parsers::{
    BismarkStats, MethrixAnnotationBySampleRow, MethrixCoverageRow, QualimapStats, SeqkitStats,
    StarStats, TrimStats,
};
use qctb::qc_summary::{
    write_excel_standard_mode, write_tsv_rnaseq, write_tsv_standard, QCSummary, QCSummaryRNA,
    ReportMode, REPORT_SCHEMA_ID,
};
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

fn sequence_stats() -> SeqkitStats {
    SeqkitStats {
        reads_raw_r1: 15,
        bases_raw_r1: 1500,
        q20_raw_r1: 91.26,
        q30_raw_r1: 81.24,
        min_len_raw_r1: 90,
        avg_len_raw_r1: 100.0,
        max_len_raw_r1: 110,
        reads_raw_r2: 15,
        bases_raw_r2: 1500,
        q20_raw_r2: 92.26,
        q30_raw_r2: 82.24,
        min_len_raw_r2: 91,
        avg_len_raw_r2: 101.0,
        max_len_raw_r2: 111,
        reads_clean_r1: 12,
        bases_clean_r1: 1200,
        q20_clean_r1: 93.26,
        q30_clean_r1: 83.24,
        min_len_clean_r1: 89,
        avg_len_clean_r1: 99.0,
        max_len_clean_r1: 109,
        reads_clean_r2: 12,
        bases_clean_r2: 1200,
        q20_clean_r2: 94.26,
        q30_clean_r2: 84.24,
        min_len_clean_r2: 88,
        avg_len_clean_r2: 98.0,
        max_len_clean_r2: 108,
        reads_raw: 30,
        bases_raw: 3000,
        reads_clean: 24,
        bases_clean: 2400,
        clean_data_ratio: 0.8,
    }
}

fn trim_stats() -> TrimStats {
    TrimStats {
        reads_with_adapter_r1: "1 (1.0%)".to_string(),
        reads_write_r1: "12".to_string(),
        bp_qc_remove_r1: "100 bp".to_string(),
        bp_write_r1: "1200 bp".to_string(),
        reads_with_adapter_r2: "1 (1.0%)".to_string(),
        reads_write_r2: "12".to_string(),
        bp_qc_remove_r2: "100 bp".to_string(),
        bp_write_r2: "1200 bp".to_string(),
    }
}

fn standard_summary(sample_id: &str, include_methrix: bool) -> QCSummary {
    let methrix_coverage = include_methrix.then_some(MethrixCoverageRow {
        sample: sample_id.to_string(),
        total_cpgs: 1000,
        covered_cpgs: 800,
        cov_1x: 800,
        cov_2x: 700,
        cov_3x: 600,
        cov_4x: 500,
        cov_5x: 400,
        cov_10x: 100,
    });
    let methrix_annotation = include_methrix.then(|| {
        let metrics = [
            ("Promoter_count", 200.0),
            ("Promoter_percent", 25.0),
            ("Exon_count", 200.0),
            ("Exon_percent", 25.0),
            ("Intron_count", 200.0),
            ("Intron_percent", 25.0),
            ("Intergenic_count", 200.0),
            ("Intergenic_percent", 25.0),
        ]
        .into_iter()
        .map(|(name, value)| (name.to_string(), value))
        .collect::<BTreeMap<_, _>>();
        MethrixAnnotationBySampleRow {
            sample: sample_id.to_string(),
            covered_cpgs: 800,
            metrics,
        }
    });

    QCSummary {
        sample_id: sample_id.to_string(),
        seqkit_stats: sequence_stats(),
        trim_stats: trim_stats(),
        bismark_stats: Some(BismarkStats {
            mapping_ratio: "80.00%".to_string(),
            total_reads_pairs: "100".to_string(),
            aligned_reads_pairs: "80".to_string(),
            aligned_reads_pairs_ratio: 0.8,
        }),
        qualimap_stats: Some(QualimapStats {
            mapping_quality: "42".to_string(),
            duplicated_reads: "10".to_string(),
            duplication_ratio: "10.00%".to_string(),
        }),
        methrix_coverage,
        methrix_annotation,
    }
}

fn rnaseq_summary() -> QCSummaryRNA {
    QCSummaryRNA {
        sample_id: "rna_sample".to_string(),
        seqkit_stats: sequence_stats(),
        trim_stats: trim_stats(),
        star_stats: StarStats {
            mapping_ratio: "75.00%".to_string(),
            total_reads: "200".to_string(),
            uniquely_mapped_reads: "150".to_string(),
            uniquely_mapped_ratio: 0.75,
        },
    }
}

fn assert_standard_golden(mode: ReportMode, summary: QCSummary, expected: &str) -> Result<()> {
    let temp_dir = TempDir::new()?;
    let output = temp_dir.path().join("report.tsv");
    write_tsv_standard(&[summary], &output, mode)?;
    assert_eq!(std::fs::read_to_string(output)?, expected);
    Ok(())
}

#[test]
fn rrbs_report_matches_golden() -> Result<()> {
    assert_standard_golden(
        ReportMode::Rrbs,
        standard_summary("rrbs_sample", true),
        include_str!("golden/rrbs.tsv.golden"),
    )
}

#[test]
fn wgbs_report_matches_golden() -> Result<()> {
    assert_standard_golden(
        ReportMode::Wgbs,
        standard_summary("wgbs_sample", false),
        include_str!("golden/wgbs.tsv.golden"),
    )
}

#[test]
fn pdx_report_matches_golden() -> Result<()> {
    assert_standard_golden(
        ReportMode::Pdx,
        standard_summary("pdx_sample", true),
        include_str!("golden/pdx.tsv.golden"),
    )
}

#[test]
fn rnaseq_report_matches_golden() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let output = temp_dir.path().join("report.tsv");
    write_tsv_rnaseq(&[rnaseq_summary()], &output)?;
    assert_eq!(
        std::fs::read_to_string(output)?,
        include_str!("golden/rnaseq.tsv.golden")
    );
    Ok(())
}

#[test]
fn excel_report_embeds_schema_mode_and_column_contract() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let output = temp_dir.path().join("report.xlsx");
    write_excel_standard_mode(
        &[standard_summary("pdx_sample", true)],
        output.to_str().unwrap(),
        ReportMode::Pdx,
    )?;

    let mut workbook = open_workbook_auto(&output)?;
    let metadata = workbook.worksheet_range("qctb_metadata")?;
    assert_eq!(
        metadata.get_value((2, 1)).unwrap().to_string(),
        REPORT_SCHEMA_ID
    );
    assert_eq!(metadata.get_value((3, 1)).unwrap().to_string(), "PDX");
    let report = workbook.worksheet_range("Report")?;
    assert_eq!(report.width(), ReportMode::Pdx.columns().len());
    assert_eq!(report.get_value((0, 0)).unwrap().to_string(), "Sample ID");
    Ok(())
}

#[test]
fn cli_pdx_smoke_handles_paths_with_spaces_and_embeds_schema() -> Result<()> {
    let temp_dir = TempDir::new()?;
    let workspace = temp_dir.path().join("analysis with spaces");
    let bismark_dir = workspace.join("bismark reports");
    std::fs::create_dir_all(&bismark_dir)?;

    let testdata = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("testdata");
    std::fs::copy(
        testdata
            .join("bismark")
            .join("0531LPHSC70203_val_1_bismark_bt2_PE_report.txt"),
        bismark_dir.join("F6703_372760_val_1_bismark_bt2_PE_report.txt"),
    )?;

    let config_path = workspace.join("config file.yaml");
    let output_path = workspace.join("qc summary.tsv");
    let config = format!(
        "SIDs:\n  - F6703_372760\nqcDir: '{}'\nqcDir_before: '{}'\nqcDir_after: '{}'\ntrimDir: '{}'\nbsmapDir: '{}'\noutDir_mCall: ''\nworkflow:\n  mode: PDX\n  species:\n    graft: human\n",
        testdata.display(),
        testdata.join("fastqcx_raw").display(),
        testdata.join("fastqcx_clean").display(),
        testdata.join("trim").display(),
        bismark_dir.display()
    );
    std::fs::write(&config_path, config)?;

    let command_output = Command::new(env!("CARGO_BIN_EXE_qctb"))
        .args([
            "--config",
            config_path.to_str().unwrap(),
            "--output",
            output_path.to_str().unwrap(),
            "--format",
            "tsv",
        ])
        .output()?;
    assert!(
        command_output.status.success(),
        "qctb failed: {}",
        String::from_utf8_lossy(&command_output.stderr)
    );
    let report = std::fs::read_to_string(output_path)?;
    assert!(report.starts_with("# qctb_schema=qctb.report/1.0.0\n# qctb_mode=PDX\n"));
    assert!(report.contains("F6703_372760"));
    Ok(())
}
