use anyhow::{bail, Context, Result};
use regex::Regex;
use std::fs::read_to_string;

#[derive(Debug, Clone)]
pub struct QualimapStats {
    pub mapping_quality: String,
    pub duplicated_reads: String,
    pub duplication_ratio: String,
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
            "Missing required Qualimap field '{}' in {}",
            field,
            file_path
        ),
        _ => bail!(
            "Duplicate Qualimap field '{}' appears {} times in {}",
            field,
            matches.len(),
            file_path
        ),
    }
}

fn parse_finite(value: &str, field: &str, file_path: &str) -> Result<f64> {
    let parsed: f64 = value
        .parse()
        .with_context(|| format!("Invalid Qualimap {} in {}: '{}'", field, file_path, value))?;
    if !parsed.is_finite() {
        bail!("Qualimap {} is not finite in {}", field, file_path);
    }
    Ok(parsed)
}

pub fn parse_qualimap_report(file_path: &str) -> Result<QualimapStats> {
    let content = read_to_string(file_path)
        .with_context(|| format!("Failed to read qualimap report: {}", file_path))?;
    let mapping_quality_regex = Regex::new(r"(?m)^\s*mean mapping quality\s*=\s*([\d.]+)\s*$")?;
    let duplicated_regex =
        Regex::new(r"(?m)^\s*number of duplicated reads(?:\s*\(estimated\))?\s*=\s*([\d,]+)\s*$")?;
    let duplication_regex = Regex::new(r"(?m)^\s*duplication rate\s*=\s*([\d.]+%?)\s*$")?;

    let mapping_quality = extract_unique(
        &content,
        &mapping_quality_regex,
        "mean mapping quality",
        file_path,
    )?;
    let mapping_quality_value = parse_finite(&mapping_quality, "mean mapping quality", file_path)?;
    if !(0.0..=255.0).contains(&mapping_quality_value) {
        bail!(
            "Qualimap mean mapping quality must be between 0 and 255 in {}, got '{}'",
            file_path,
            mapping_quality
        );
    }

    let duplicated_reads = extract_unique(
        &content,
        &duplicated_regex,
        "number of duplicated reads",
        file_path,
    )?;
    duplicated_reads.parse::<u64>().with_context(|| {
        format!(
            "Invalid Qualimap duplicated-read count in {}: '{}'",
            file_path, duplicated_reads
        )
    })?;

    let duplication_ratio =
        extract_unique(&content, &duplication_regex, "duplication rate", file_path)?;
    let (numeric_ratio, is_percent) = duplication_ratio
        .strip_suffix('%')
        .map_or((duplication_ratio.as_str(), false), |value| (value, true));
    let ratio_value = parse_finite(numeric_ratio, "duplication rate", file_path)?;
    let valid_range = if is_percent { 0.0..=100.0 } else { 0.0..=1.0 };
    if !valid_range.contains(&ratio_value) {
        bail!(
            "Qualimap duplication rate has invalid range in {}: '{}'",
            file_path,
            duplication_ratio
        );
    }

    Ok(QualimapStats {
        mapping_quality,
        duplicated_reads,
        duplication_ratio,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::PathBuf;
    use tempfile::NamedTempFile;

    fn write_valid_report(file: &mut NamedTempFile) -> Result<()> {
        writeln!(file, "qualimap report")?;
        writeln!(file, "mean mapping quality = 60")?;
        writeln!(file, "number of duplicated reads = 50000")?;
        writeln!(file, "duplication rate = 0.0526")?;
        Ok(())
    }

    #[test]
    fn test_parse_qualimap_report() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        write_valid_report(&mut temp_file)?;
        let stats = parse_qualimap_report(temp_file.path().to_str().unwrap())?;
        assert_eq!(stats.mapping_quality, "60");
        assert_eq!(stats.duplicated_reads, "50000");
        assert_eq!(stats.duplication_ratio, "0.0526");
        Ok(())
    }

    #[test]
    fn real_report_preserves_percent_and_comma_counts() -> Result<()> {
        let file = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("testdata/qualimap/F6703_372760_human/genome_results.txt");
        let stats = parse_qualimap_report(file.to_str().unwrap())?;
        assert_eq!(stats.mapping_quality, "15.6988");
        assert_eq!(stats.duplicated_reads, "35998499");
        assert_eq!(stats.duplication_ratio, "51.22%");
        Ok(())
    }

    #[test]
    fn test_missing_required_field_fails() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, "mean mapping quality = 60")?;
        writeln!(temp_file, "number of duplicated reads = 50000")?;
        assert!(parse_qualimap_report(temp_file.path().to_str().unwrap()).is_err());
        Ok(())
    }

    #[test]
    fn duplicate_and_out_of_range_fields_are_rejected() -> Result<()> {
        let mut duplicate_file = NamedTempFile::new()?;
        write_valid_report(&mut duplicate_file)?;
        writeln!(duplicate_file, "mean mapping quality = 60")?;
        assert!(parse_qualimap_report(duplicate_file.path().to_str().unwrap()).is_err());

        let mut invalid_rate_file = NamedTempFile::new()?;
        writeln!(invalid_rate_file, "mean mapping quality = 60")?;
        writeln!(invalid_rate_file, "number of duplicated reads = 50000")?;
        writeln!(invalid_rate_file, "duplication rate = 1.5")?;
        assert!(parse_qualimap_report(invalid_rate_file.path().to_str().unwrap()).is_err());
        Ok(())
    }
}
