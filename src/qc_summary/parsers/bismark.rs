use anyhow::{bail, Context, Result};
use regex::Regex;
use std::fs::read_to_string;

#[derive(Debug, Clone)]
pub struct BismarkStats {
    pub mapping_ratio: String,
    pub total_reads_pairs: String,
    pub aligned_reads_pairs: String,
    pub aligned_reads_pairs_ratio: f64,
}

fn extract_unique(content: &str, regex: &Regex, field: &str, file_path: &str) -> Result<String> {
    let matches = regex
        .captures_iter(content)
        .filter_map(|captures| captures.get(1))
        .map(|value| value.as_str().replace(',', ""))
        .collect::<Vec<_>>();
    match matches.as_slice() {
        [value] => Ok(value.clone()),
        [] => bail!(
            "Missing required Bismark field '{}' in {}",
            field,
            file_path
        ),
        _ => bail!(
            "Duplicate Bismark field '{}' appears {} times in {}",
            field,
            matches.len(),
            file_path
        ),
    }
}

fn parse_percent(value: &str, field: &str, file_path: &str) -> Result<f64> {
    let parsed: f64 = value
        .parse()
        .with_context(|| format!("Invalid Bismark {} in {}: '{}'", field, file_path, value))?;
    if !parsed.is_finite() || !(0.0..=100.0).contains(&parsed) {
        bail!(
            "Bismark {} must be a finite percentage from 0 to 100 in {}, got '{}'",
            field,
            file_path,
            value
        );
    }
    Ok(parsed)
}

pub fn parse_bismark_report(file_path: &str) -> Result<BismarkStats> {
    let content = read_to_string(file_path)
        .with_context(|| format!("Failed to read bismark report: {}", file_path))?;
    let mapping_regex = Regex::new(r"(?m)^\s*Mapping efficiency:\s*([\d.]+)%\s*$")?;
    let total_regex = Regex::new(r"(?m)^\s*Sequence pairs analysed in total:\s+([\d,]+)\s*$")?;
    let aligned_regex = Regex::new(
        r"(?m)^\s*Number of paired-end alignments with a unique best hit:\s+([\d,]+)\s*$",
    )?;

    let mapping_value = extract_unique(&content, &mapping_regex, "Mapping efficiency", file_path)?;
    let mapping_percent = parse_percent(&mapping_value, "mapping efficiency", file_path)?;
    let total_reads_pairs = extract_unique(
        &content,
        &total_regex,
        "Sequence pairs analysed in total",
        file_path,
    )?;
    let aligned_reads_pairs = extract_unique(
        &content,
        &aligned_regex,
        "Number of paired-end alignments with a unique best hit",
        file_path,
    )?;

    let total_count: u64 = total_reads_pairs
        .parse()
        .with_context(|| format!("Invalid total read-pair count in {}", file_path))?;
    let aligned_count: u64 = aligned_reads_pairs
        .parse()
        .with_context(|| format!("Invalid aligned read-pair count in {}", file_path))?;
    if total_count == 0 {
        bail!("Bismark total read-pair count is zero in {}", file_path);
    }
    if aligned_count > total_count {
        bail!(
            "Bismark aligned read-pair count {} exceeds total {} in {}",
            aligned_count,
            total_count,
            file_path
        );
    }
    let aligned_reads_pairs_ratio = aligned_count as f64 / total_count as f64;
    let calculated_percent = aligned_reads_pairs_ratio * 100.0;
    if (calculated_percent - mapping_percent).abs() > 0.05 {
        bail!(
            "Bismark mapping efficiency {}% disagrees with aligned/total calculation {:.4}% in {}",
            mapping_percent,
            calculated_percent,
            file_path
        );
    }

    Ok(BismarkStats {
        mapping_ratio: format!("{}%", mapping_value),
        total_reads_pairs,
        aligned_reads_pairs,
        aligned_reads_pairs_ratio,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_valid_report(file: &mut NamedTempFile) -> Result<()> {
        writeln!(file, "Bismark PE report")?;
        writeln!(file, "Mapping efficiency: 95.23%")?;
        writeln!(file, "Sequence pairs analysed in total: 1000000")?;
        writeln!(
            file,
            "Number of paired-end alignments with a unique best hit: 952300"
        )?;
        Ok(())
    }

    #[test]
    fn test_parse_bismark_report() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        write_valid_report(&mut temp_file)?;
        let stats = parse_bismark_report(temp_file.path().to_str().unwrap())?;
        assert_eq!(stats.mapping_ratio, "95.23%");
        assert_eq!(stats.total_reads_pairs, "1000000");
        assert_eq!(stats.aligned_reads_pairs, "952300");
        assert!((stats.aligned_reads_pairs_ratio - 0.9523).abs() < 0.0001);
        Ok(())
    }

    #[test]
    fn test_missing_required_field_fails() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, "Mapping efficiency: 95.23%")?;
        writeln!(temp_file, "Sequence pairs analysed in total: 1000000")?;
        assert!(parse_bismark_report(temp_file.path().to_str().unwrap()).is_err());
        Ok(())
    }

    #[test]
    fn duplicate_fields_are_rejected() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        write_valid_report(&mut temp_file)?;
        writeln!(temp_file, "Mapping efficiency: 95.23%")?;
        assert!(parse_bismark_report(temp_file.path().to_str().unwrap()).is_err());
        Ok(())
    }

    #[test]
    fn contradictory_counts_and_percentage_are_rejected() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, "Mapping efficiency: 95.23%")?;
        writeln!(temp_file, "Sequence pairs analysed in total: 100")?;
        writeln!(
            temp_file,
            "Number of paired-end alignments with a unique best hit: 101"
        )?;
        assert!(parse_bismark_report(temp_file.path().to_str().unwrap()).is_err());

        let mut inconsistent_file = NamedTempFile::new()?;
        writeln!(inconsistent_file, "Mapping efficiency: 50.0%")?;
        writeln!(inconsistent_file, "Sequence pairs analysed in total: 100")?;
        writeln!(
            inconsistent_file,
            "Number of paired-end alignments with a unique best hit: 90"
        )?;
        assert!(parse_bismark_report(inconsistent_file.path().to_str().unwrap()).is_err());
        Ok(())
    }
}
