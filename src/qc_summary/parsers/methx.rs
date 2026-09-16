use crate::qc_summary::schema::CONTRACT_ANNOTATION_METRICS;
use anyhow::{bail, Context, Result};
use calamine::{open_workbook_auto, Data, Reader};
use std::collections::{BTreeMap, HashSet};

#[derive(Debug, Clone, PartialEq)]
pub struct MethrixCoverageRow {
    pub sample: String,
    pub total_cpgs: u64,
    pub covered_cpgs: u64,
    pub cov_1x: u64,
    pub cov_2x: u64,
    pub cov_3x: u64,
    pub cov_4x: u64,
    pub cov_5x: u64,
    pub cov_10x: u64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct MethrixAnnotationBySampleRow {
    pub sample: String,
    pub covered_cpgs: u64,
    pub metrics: BTreeMap<String, f64>,
}

fn cell_as_string(cell: Option<&Data>) -> String {
    match cell {
        Some(value) => value.to_string().trim().to_string(),
        None => String::new(),
    }
}

fn cell_as_u64(cell: Option<&Data>, field: &str) -> Result<u64> {
    let raw = cell_as_string(cell);
    let value: f64 = raw
        .parse()
        .with_context(|| format!("Failed to parse '{}' as number: '{}'", field, raw))?;
    if !value.is_finite() || value < 0.0 || value.fract() != 0.0 || value > u64::MAX as f64 {
        bail!(
            "Field '{}' must be a finite non-negative integer, got '{}'",
            field,
            raw
        );
    }
    Ok(value as u64)
}

fn cell_as_f64(cell: Option<&Data>, field: &str) -> Result<f64> {
    let raw = cell_as_string(cell);
    let value: f64 = raw
        .parse()
        .with_context(|| format!("Failed to parse '{}' as number: '{}'", field, raw))?;
    if !value.is_finite() || value < 0.0 {
        bail!(
            "Field '{}' must be finite and non-negative, got '{}'",
            field,
            raw
        );
    }
    Ok(value)
}

fn validated_headers(header: &[Data], sheet: &str) -> Result<Vec<String>> {
    let headers = header
        .iter()
        .map(|cell| cell.to_string().trim().to_string())
        .collect::<Vec<_>>();
    let mut seen = HashSet::new();
    for name in &headers {
        if name.is_empty() {
            bail!("{} contains an empty column name", sheet);
        }
        if !seen.insert(name) {
            bail!("{} contains duplicate column '{}'", sheet, name);
        }
    }
    Ok(headers)
}

fn find_column(headers: &[String], name: &str, sheet: &str) -> Result<usize> {
    headers
        .iter()
        .position(|header| header == name)
        .with_context(|| format!("Missing required column '{}' in {}", name, sheet))
}

fn validate_coverage(row: &MethrixCoverageRow, path: &str) -> Result<()> {
    if row.covered_cpgs > row.total_cpgs {
        bail!(
            "Methrix covered CpGs exceed total CpGs for sample '{}' in {}",
            row.sample,
            path
        );
    }
    let thresholds = [
        row.cov_1x,
        row.cov_2x,
        row.cov_3x,
        row.cov_4x,
        row.cov_5x,
        row.cov_10x,
    ];
    if thresholds[0] > row.covered_cpgs || thresholds.windows(2).any(|pair| pair[1] > pair[0]) {
        bail!(
            "Methrix coverage thresholds are not monotonically decreasing for sample '{}' in {}",
            row.sample,
            path
        );
    }
    Ok(())
}

pub fn parse_methx_coverage_xlsx(path: &str) -> Result<Vec<MethrixCoverageRow>> {
    let mut workbook =
        open_workbook_auto(path).with_context(|| format!("Failed to open XLSX: {}", path))?;
    let range = workbook
        .worksheet_range("Sheet1")
        .with_context(|| format!("Failed to read Sheet1 in {}", path))?;

    let mut rows = range.rows();
    let header = rows.next().context("Coverage sheet is empty")?;
    let headers = validated_headers(header, "Sheet1")?;
    let sample_index = find_column(&headers, "Sample", "Sheet1")?;
    let total_index = find_column(&headers, "Total CpGs", "Sheet1")?;
    let covered_index = find_column(&headers, "Covered CpGs", "Sheet1")?;
    let coverage_indices = ["1X", "2X", "3X", "4X", "5X", "10X"]
        .map(|name| find_column(&headers, name, "Sheet1"))
        .into_iter()
        .collect::<Result<Vec<_>>>()?;

    let mut output = Vec::new();
    let mut seen_samples = HashSet::new();
    for row in rows {
        let sample = cell_as_string(row.get(sample_index));
        if sample.is_empty() {
            continue;
        }
        if !seen_samples.insert(sample.clone()) {
            bail!("Duplicate Methrix coverage sample '{}' in {}", sample, path);
        }
        let parsed = MethrixCoverageRow {
            sample,
            total_cpgs: cell_as_u64(row.get(total_index), "Total CpGs")?,
            covered_cpgs: cell_as_u64(row.get(covered_index), "Covered CpGs")?,
            cov_1x: cell_as_u64(row.get(coverage_indices[0]), "1X")?,
            cov_2x: cell_as_u64(row.get(coverage_indices[1]), "2X")?,
            cov_3x: cell_as_u64(row.get(coverage_indices[2]), "3X")?,
            cov_4x: cell_as_u64(row.get(coverage_indices[3]), "4X")?,
            cov_5x: cell_as_u64(row.get(coverage_indices[4]), "5X")?,
            cov_10x: cell_as_u64(row.get(coverage_indices[5]), "10X")?,
        };
        validate_coverage(&parsed, path)?;
        output.push(parsed);
    }
    Ok(output)
}

pub fn parse_methx_annotation_by_sample_xlsx(
    path: &str,
) -> Result<Vec<MethrixAnnotationBySampleRow>> {
    let mut workbook =
        open_workbook_auto(path).with_context(|| format!("Failed to open XLSX: {}", path))?;
    let sheet_name = "ChIPseeker_By_Sample";
    let range = workbook
        .worksheet_range(sheet_name)
        .with_context(|| format!("Failed to read {} in {}", sheet_name, path))?;

    let mut rows = range.rows();
    let header = rows.next().context("Annotation-by-sample sheet is empty")?;
    let headers = validated_headers(header, sheet_name)?;
    let sample_index = find_column(&headers, "sample", sheet_name)?;
    let covered_index = find_column(&headers, "covered_cpgs", sheet_name)?;
    for metric in CONTRACT_ANNOTATION_METRICS {
        find_column(&headers, metric, sheet_name)?;
    }

    let metric_indices = headers
        .iter()
        .enumerate()
        .filter(|(index, _)| *index != sample_index && *index != covered_index)
        .map(|(index, header)| (index, header.clone()))
        .collect::<Vec<_>>();

    let mut output = Vec::new();
    let mut seen_samples = HashSet::new();
    for row in rows {
        let sample = cell_as_string(row.get(sample_index));
        if sample.is_empty() {
            continue;
        }
        if !seen_samples.insert(sample.clone()) {
            bail!(
                "Duplicate Methrix annotation sample '{}' in {}",
                sample,
                path
            );
        }
        let mut metrics = BTreeMap::new();
        for (index, name) in &metric_indices {
            let value = cell_as_f64(row.get(*index), name)?;
            if name.ends_with("_count") && value.fract() != 0.0 {
                bail!(
                    "Methrix count metric '{}' must be an integer in {}",
                    name,
                    path
                );
            }
            if name.ends_with("_percent") && value > 100.0 {
                bail!("Methrix percent metric '{}' exceeds 100 in {}", name, path);
            }
            metrics.insert(name.clone(), value);
        }
        output.push(MethrixAnnotationBySampleRow {
            sample,
            covered_cpgs: cell_as_u64(row.get(covered_index), "covered_cpgs")?,
            metrics,
        });
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_xlsxwriter::Workbook;
    use tempfile::TempDir;

    fn write_coverage(path: &std::path::Path, invalid_order: bool) -> Result<()> {
        let mut workbook = Workbook::new();
        let sheet = workbook.add_worksheet();
        for (column, header) in [
            "Sample",
            "Total CpGs",
            "Covered CpGs",
            "1X",
            "2X",
            "3X",
            "4X",
            "5X",
            "10X",
        ]
        .iter()
        .enumerate()
        {
            sheet.write_string(0, column as u16, *header)?;
        }
        sheet.write_string(1, 0, "sample1_nsort.bismark.cov")?;
        for (column, value) in [100.0, 80.0, 80.0, 70.0, 60.0, 50.0, 40.0, 10.0]
            .iter()
            .enumerate()
        {
            let value = if invalid_order && column == 3 {
                90.0
            } else {
                *value
            };
            sheet.write_number(1, (column + 1) as u16, value)?;
        }
        workbook.save(path)?;
        Ok(())
    }

    fn write_annotation(path: &std::path::Path, omit_last_metric: bool) -> Result<()> {
        let mut workbook = Workbook::new();
        let sheet = workbook.add_worksheet();
        sheet.set_name("ChIPseeker_By_Sample")?;
        sheet.write_string(0, 0, "sample")?;
        sheet.write_string(0, 1, "covered_cpgs")?;
        let metrics: &[&str] = if omit_last_metric {
            &CONTRACT_ANNOTATION_METRICS[..7]
        } else {
            &CONTRACT_ANNOTATION_METRICS
        };
        for (column, metric) in metrics.iter().enumerate() {
            sheet.write_string(0, (column + 2) as u16, *metric)?;
        }
        sheet.write_string(1, 0, "sample1_nsort.bismark.cov")?;
        sheet.write_number(1, 1, 80)?;
        for (column, metric) in metrics.iter().enumerate() {
            let value = if metric.ends_with("_count") {
                10.0
            } else {
                12.5
            };
            sheet.write_number(1, (column + 2) as u16, value)?;
        }
        workbook.save(path)?;
        Ok(())
    }

    #[test]
    fn integer_cells_reject_fractional_and_negative_values() {
        assert!(cell_as_u64(Some(&Data::Float(1.5)), "count").is_err());
        assert!(cell_as_u64(Some(&Data::Float(-1.0)), "count").is_err());
        assert_eq!(cell_as_u64(Some(&Data::Float(2.0)), "count").unwrap(), 2);
    }

    #[test]
    fn parses_runtime_generated_coverage_report() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let path = temp_dir.path().join("coverage.xlsx");
        write_coverage(&path, false)?;
        let rows = parse_methx_coverage_xlsx(path.to_str().unwrap())?;
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].sample, "sample1_nsort.bismark.cov");
        assert_eq!(rows[0].covered_cpgs, 80);
        assert_eq!(rows[0].cov_10x, 10);
        Ok(())
    }

    #[test]
    fn inconsistent_coverage_thresholds_fail() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let path = temp_dir.path().join("coverage.xlsx");
        write_coverage(&path, true)?;
        assert!(parse_methx_coverage_xlsx(path.to_str().unwrap()).is_err());
        Ok(())
    }

    #[test]
    fn parses_runtime_generated_annotation_report() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let path = temp_dir.path().join("annotation.xlsx");
        write_annotation(&path, false)?;
        let rows = parse_methx_annotation_by_sample_xlsx(path.to_str().unwrap())?;
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].covered_cpgs, 80);
        assert_eq!(rows[0].metrics.get("Promoter_count"), Some(&10.0));
        Ok(())
    }

    #[test]
    fn missing_contracted_annotation_metric_fails() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let path = temp_dir.path().join("annotation.xlsx");
        write_annotation(&path, true)?;
        assert!(parse_methx_annotation_by_sample_xlsx(path.to_str().unwrap()).is_err());
        Ok(())
    }
}
