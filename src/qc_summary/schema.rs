use crate::qc_summary::{QCSummary, QCSummaryRNA};
use anyhow::{bail, Context, Result};

pub const REPORT_SCHEMA_NAME: &str = "qctb.report";
pub const REPORT_SCHEMA_VERSION: &str = "1.0.0";
pub const REPORT_SCHEMA_ID: &str = "qctb.report/1.0.0";
pub const REPORT_MISSING_VALUE: &str = "N/A";

pub const CONTRACT_ANNOTATION_METRICS: [&str; 8] = [
    "Promoter_count",
    "Promoter_percent",
    "Exon_count",
    "Exon_percent",
    "Intron_count",
    "Intron_percent",
    "Intergenic_count",
    "Intergenic_percent",
];

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReportMode {
    Rrbs,
    Wgbs,
    RnaSeq,
    Pdx,
    Standard,
}

impl ReportMode {
    pub fn from_workflow_mode(mode: Option<&str>, rnaseq: bool) -> Self {
        if rnaseq {
            return Self::RnaSeq;
        }
        match mode.map(str::trim).map(str::to_ascii_uppercase).as_deref() {
            Some("RRBS") => Self::Rrbs,
            Some("WGBS") | Some("BSSEQ") => Self::Wgbs,
            Some("PDX") | Some("BEAVERPDX") | Some("BEAVERRNASEQPDX") => Self::Pdx,
            Some("RNASEQ") | Some("RNA-SEQ") => Self::RnaSeq,
            _ => Self::Standard,
        }
    }

    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Rrbs => "RRBS",
            Self::Wgbs => "WGBS",
            Self::RnaSeq => "RNA-seq",
            Self::Pdx => "PDX",
            Self::Standard => "standard",
        }
    }

    pub const fn columns(self) -> &'static [ReportColumnSpec] {
        match self {
            Self::RnaSeq => RNA_COLUMNS,
            Self::Rrbs | Self::Wgbs | Self::Pdx | Self::Standard => STANDARD_COLUMNS,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReportColumnType {
    Text,
    Integer,
    Decimal,
}

impl ReportColumnType {
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Text => "text",
            Self::Integer => "integer",
            Self::Decimal => "decimal",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReportColumnSpec {
    pub key: &'static str,
    pub excel_header: &'static str,
    pub column_type: ReportColumnType,
    pub decimal_places: Option<u8>,
}

const fn text(key: &'static str, excel_header: &'static str) -> ReportColumnSpec {
    ReportColumnSpec {
        key,
        excel_header,
        column_type: ReportColumnType::Text,
        decimal_places: None,
    }
}

const fn integer(key: &'static str, excel_header: &'static str) -> ReportColumnSpec {
    ReportColumnSpec {
        key,
        excel_header,
        column_type: ReportColumnType::Integer,
        decimal_places: None,
    }
}

const fn decimal(
    key: &'static str,
    excel_header: &'static str,
    decimal_places: u8,
) -> ReportColumnSpec {
    ReportColumnSpec {
        key,
        excel_header,
        column_type: ReportColumnType::Decimal,
        decimal_places: Some(decimal_places),
    }
}

pub const STANDARD_COLUMNS: &[ReportColumnSpec] = &[
    text("sample_id", "Sample ID"),
    integer("reads_raw", "Reads Raw"),
    integer("bases_raw", "Bases Raw"),
    integer("reads_clean", "Reads Clean"),
    integer("bases_clean", "Bases Clean"),
    decimal("clean_data_ratio", "Clean Data Ratio", 4),
    decimal("q20_raw_r1", "Q20 Raw R1 (%)", 1),
    decimal("q30_raw_r1", "Q30 Raw R1 (%)", 1),
    decimal("avg_len_raw_r1", "Avg Len Raw R1", 1),
    decimal("q20_raw_r2", "Q20 Raw R2 (%)", 1),
    decimal("q30_raw_r2", "Q30 Raw R2 (%)", 1),
    decimal("avg_len_raw_r2", "Avg Len Raw R2", 1),
    decimal("q20_clean_r1", "Q20 Clean R1 (%)", 1),
    decimal("q30_clean_r1", "Q30 Clean R1 (%)", 1),
    decimal("avg_len_clean_r1", "Avg Len Clean R1", 1),
    decimal("q20_clean_r2", "Q20 Clean R2 (%)", 1),
    decimal("q30_clean_r2", "Q30 Clean R2 (%)", 1),
    decimal("avg_len_clean_r2", "Avg Len Clean R2", 1),
    text("mapping_ratio", "Mapping Ratio (%)"),
    text("total_reads_pairs", "Total Read Pairs"),
    text("aligned_reads_pairs", "Aligned Read Pairs"),
    decimal("aligned_ratio", "Aligned Ratio", 4),
    text("mapping_quality", "Mapping Quality"),
    text("duplicated_reads", "Duplicated Reads"),
    text("duplication_ratio", "Duplication Rate"),
    integer("methrix_total_cpgs", "Methrix Total CpGs"),
    integer("methrix_covered_cpgs", "Methrix Covered CpGs"),
    integer("methrix_1x", "Methrix 1X"),
    integer("methrix_2x", "Methrix 2X"),
    integer("methrix_3x", "Methrix 3X"),
    integer("methrix_4x", "Methrix 4X"),
    integer("methrix_5x", "Methrix 5X"),
    integer("methrix_10x", "Methrix 10X"),
    integer("methrix_ann_covered_cpgs", "Methrix Ann Covered CpGs"),
    decimal("methrix_promoter_count", "Methrix Promoter Count", 0),
    decimal("methrix_promoter_percent", "Methrix Promoter Percent", 6),
    decimal("methrix_exon_count", "Methrix Exon Count", 0),
    decimal("methrix_exon_percent", "Methrix Exon Percent", 6),
    decimal("methrix_intron_count", "Methrix Intron Count", 0),
    decimal("methrix_intron_percent", "Methrix Intron Percent", 6),
    decimal("methrix_intergenic_count", "Methrix Intergenic Count", 0),
    decimal(
        "methrix_intergenic_percent",
        "Methrix Intergenic Percent",
        6,
    ),
];

pub const RNA_COLUMNS: &[ReportColumnSpec] = &[
    text("sample_id", "Sample ID"),
    integer("reads_raw", "Reads Raw"),
    integer("bases_raw", "Bases Raw"),
    integer("reads_clean", "Reads Clean"),
    integer("bases_clean", "Bases Clean"),
    decimal("clean_data_ratio", "Clean Data Ratio", 4),
    decimal("q20_raw_r1", "Q20 Raw R1 (%)", 1),
    decimal("q30_raw_r1", "Q30 Raw R1 (%)", 1),
    decimal("avg_len_raw_r1", "Avg Len Raw R1", 1),
    decimal("q20_raw_r2", "Q20 Raw R2 (%)", 1),
    decimal("q30_raw_r2", "Q30 Raw R2 (%)", 1),
    decimal("avg_len_raw_r2", "Avg Len Raw R2", 1),
    decimal("q20_clean_r1", "Q20 Clean R1 (%)", 1),
    decimal("q30_clean_r1", "Q30 Clean R1 (%)", 1),
    decimal("avg_len_clean_r1", "Avg Len Clean R1", 1),
    decimal("q20_clean_r2", "Q20 Clean R2 (%)", 1),
    decimal("q30_clean_r2", "Q30 Clean R2 (%)", 1),
    decimal("avg_len_clean_r2", "Avg Len Clean R2", 1),
    text("mapping_ratio", "Mapping Ratio (%)"),
    text("total_reads", "Total Reads"),
    text("uniquely_mapped_reads", "Uniquely Mapped Reads"),
    decimal("uniquely_mapped_ratio", "Uniquely Mapped Ratio", 4),
];

#[derive(Clone, Debug, PartialEq)]
pub enum ReportCell {
    Text(String),
    Integer(u64),
    Decimal(f64),
    Missing,
}

impl ReportCell {
    pub fn tsv_value(&self, specification: &ReportColumnSpec) -> Result<String> {
        match (self, specification.column_type) {
            (Self::Text(value), ReportColumnType::Text) => {
                if value.contains(['\t', '\r', '\n']) {
                    bail!("TSV text value contains a tab or newline: {value:?}");
                }
                Ok(value.clone())
            }
            (Self::Integer(value), ReportColumnType::Integer) => Ok(value.to_string()),
            (Self::Decimal(value), ReportColumnType::Decimal) => {
                if !value.is_finite() {
                    bail!(
                        "Report decimal is not finite for column '{}'",
                        specification.key
                    );
                }
                let places = specification.decimal_places.unwrap_or(0) as usize;
                Ok(format!("{value:.places$}"))
            }
            (Self::Missing, _) => Ok(REPORT_MISSING_VALUE.to_string()),
            _ => bail!(
                "Report value type does not match schema column '{}' ({})",
                specification.key,
                specification.column_type.as_str()
            ),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct ReportRow {
    pub cells: Vec<ReportCell>,
}

impl ReportRow {
    pub fn validate(&self, mode: ReportMode) -> Result<()> {
        let columns = mode.columns();
        if self.cells.len() != columns.len() {
            bail!(
                "Report row has {} cells but {} mode requires {}",
                self.cells.len(),
                mode.as_str(),
                columns.len()
            );
        }
        for (cell, column) in self.cells.iter().zip(columns) {
            cell.tsv_value(column)?;
        }
        Ok(())
    }
}

fn common_cells(
    sample_id: &str,
    stats: &crate::qc_summary::parsers::SeqkitStats,
) -> Vec<ReportCell> {
    vec![
        ReportCell::Text(sample_id.to_string()),
        ReportCell::Integer(stats.reads_raw),
        ReportCell::Integer(stats.bases_raw),
        ReportCell::Integer(stats.reads_clean),
        ReportCell::Integer(stats.bases_clean),
        ReportCell::Decimal(stats.clean_data_ratio),
        ReportCell::Decimal(stats.q20_raw_r1),
        ReportCell::Decimal(stats.q30_raw_r1),
        ReportCell::Decimal(stats.avg_len_raw_r1),
        ReportCell::Decimal(stats.q20_raw_r2),
        ReportCell::Decimal(stats.q30_raw_r2),
        ReportCell::Decimal(stats.avg_len_raw_r2),
        ReportCell::Decimal(stats.q20_clean_r1),
        ReportCell::Decimal(stats.q30_clean_r1),
        ReportCell::Decimal(stats.avg_len_clean_r1),
        ReportCell::Decimal(stats.q20_clean_r2),
        ReportCell::Decimal(stats.q30_clean_r2),
        ReportCell::Decimal(stats.avg_len_clean_r2),
    ]
}

pub fn standard_report_rows(summaries: &[QCSummary], mode: ReportMode) -> Result<Vec<ReportRow>> {
    if mode == ReportMode::RnaSeq {
        bail!("RNA-seq mode cannot be used with standard QC summaries");
    }
    summaries
        .iter()
        .map(|summary| {
            let mut cells = common_cells(&summary.sample_id, &summary.seqkit_stats);
            if let Some(stats) = &summary.bismark_stats {
                cells.extend([
                    ReportCell::Text(stats.mapping_ratio.clone()),
                    ReportCell::Text(stats.total_reads_pairs.clone()),
                    ReportCell::Text(stats.aligned_reads_pairs.clone()),
                    ReportCell::Decimal(stats.aligned_reads_pairs_ratio),
                ]);
            } else {
                cells.extend(std::iter::repeat_n(ReportCell::Missing, 4));
            }
            if let Some(stats) = &summary.qualimap_stats {
                cells.extend([
                    ReportCell::Text(stats.mapping_quality.clone()),
                    ReportCell::Text(stats.duplicated_reads.clone()),
                    ReportCell::Text(stats.duplication_ratio.clone()),
                ]);
            } else {
                cells.extend(std::iter::repeat_n(ReportCell::Missing, 3));
            }
            if let Some(stats) = &summary.methrix_coverage {
                cells.extend([
                    ReportCell::Integer(stats.total_cpgs),
                    ReportCell::Integer(stats.covered_cpgs),
                    ReportCell::Integer(stats.cov_1x),
                    ReportCell::Integer(stats.cov_2x),
                    ReportCell::Integer(stats.cov_3x),
                    ReportCell::Integer(stats.cov_4x),
                    ReportCell::Integer(stats.cov_5x),
                    ReportCell::Integer(stats.cov_10x),
                ]);
            } else {
                cells.extend(std::iter::repeat_n(ReportCell::Missing, 8));
            }
            if let Some(stats) = &summary.methrix_annotation {
                cells.push(ReportCell::Integer(stats.covered_cpgs));
                for metric in CONTRACT_ANNOTATION_METRICS {
                    cells.push(ReportCell::Decimal(
                        stats.metrics.get(metric).copied().with_context(|| {
                            format!(
                                "Methrix annotation for sample '{}' is missing contracted metric '{}'",
                                summary.sample_id, metric
                            )
                        })?,
                    ));
                }
            } else {
                cells.extend(std::iter::repeat_n(ReportCell::Missing, 9));
            }
            let row = ReportRow { cells };
            row.validate(mode)?;
            Ok(row)
        })
        .collect()
}

pub fn rnaseq_report_rows(summaries: &[QCSummaryRNA]) -> Result<Vec<ReportRow>> {
    summaries
        .iter()
        .map(|summary| {
            let mut cells = common_cells(&summary.sample_id, &summary.seqkit_stats);
            cells.extend([
                ReportCell::Text(summary.star_stats.mapping_ratio.clone()),
                ReportCell::Text(summary.star_stats.total_reads.clone()),
                ReportCell::Text(summary.star_stats.uniquely_mapped_reads.clone()),
                ReportCell::Decimal(summary.star_stats.uniquely_mapped_ratio),
            ]);
            let row = ReportRow { cells };
            row.validate(ReportMode::RnaSeq)?;
            Ok(row)
        })
        .collect()
}
