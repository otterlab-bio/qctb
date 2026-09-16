use anyhow::{bail, Context, Result};
use regex::Regex;
use std::fs::read_to_string;

#[derive(Debug, Clone)]
pub struct StarStats {
    pub mapping_ratio: String,
    pub total_reads: String,
    pub uniquely_mapped_reads: String,
    pub uniquely_mapped_ratio: f64,
}

fn extract_unique(content: &str, regex: &Regex, field: &str, file_path: &str) -> Result<String> {
    let matches = regex
        .captures_iter(content)
        .filter_map(|captures| captures.get(1))
        .map(|value| value.as_str().replace(',', ""))
        .collect::<Vec<_>>();
    match matches.as_slice() {
        [value] => Ok(value.clone()),
        [] => bail!("Missing required STAR field '{}' in {}", field, file_path),
        _ => bail!(
            "Duplicate STAR field '{}' appears {} times in {}",
            field,
            matches.len(),
            file_path
        ),
    }
}

fn parse_percent(value: &str, file_path: &str) -> Result<f64> {
    let parsed: f64 = value.parse().with_context(|| {
        format!(
            "Invalid STAR mapping percentage in {}: '{}'",
            file_path, value
        )
    })?;
    if !parsed.is_finite() || !(0.0..=100.0).contains(&parsed) {
        bail!(
            "STAR mapping percentage must be finite and from 0 to 100 in {}, got '{}'",
            file_path,
            value
        );
    }
    Ok(parsed)
}

pub fn parse_star_log(file_path: &str) -> Result<StarStats> {
    let content = read_to_string(file_path)
        .with_context(|| format!("Failed to read STAR log file: {}", file_path))?;
    let uniquely_percent_regex =
        Regex::new(r"(?m)^\s*\|?\s*Uniquely mapped reads %\s*\|\s*([\d.]+)%\s*\|?\s*$")?;
    let input_reads_regex =
        Regex::new(r"(?m)^\s*\|?\s*Number of input reads\s*\|\s*([\d,]+)\s*\|?\s*$")?;
    let uniquely_number_regex =
        Regex::new(r"(?m)^\s*\|?\s*Uniquely mapped reads number\s*\|\s*([\d,]+)\s*\|?\s*$")?;

    let mapping_value = extract_unique(
        &content,
        &uniquely_percent_regex,
        "Uniquely mapped reads %",
        file_path,
    )?;
    let mapping_percent = parse_percent(&mapping_value, file_path)?;
    let total_reads = extract_unique(
        &content,
        &input_reads_regex,
        "Number of input reads",
        file_path,
    )?;
    let uniquely_mapped_reads = extract_unique(
        &content,
        &uniquely_number_regex,
        "Uniquely mapped reads number",
        file_path,
    )?;

    let total_count: u64 = total_reads
        .parse()
        .with_context(|| format!("Invalid input-read count in {}", file_path))?;
    let uniquely_mapped_count: u64 = uniquely_mapped_reads
        .parse()
        .with_context(|| format!("Invalid uniquely mapped read count in {}", file_path))?;
    if total_count == 0 {
        bail!("STAR input-read count is zero in {}", file_path);
    }
    if uniquely_mapped_count > total_count {
        bail!(
            "STAR uniquely mapped read count {} exceeds total {} in {}",
            uniquely_mapped_count,
            total_count,
            file_path
        );
    }
    let uniquely_mapped_ratio = uniquely_mapped_count as f64 / total_count as f64;
    let calculated_percent = uniquely_mapped_ratio * 100.0;
    if (calculated_percent - mapping_percent).abs() > 0.05 {
        bail!(
            "STAR mapping percentage {}% disagrees with unique/total calculation {:.4}% in {}",
            mapping_percent,
            calculated_percent,
            file_path
        );
    }

    Ok(StarStats {
        mapping_ratio: format!("{}%", mapping_value),
        total_reads,
        uniquely_mapped_reads,
        uniquely_mapped_ratio,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_valid_log(file: &mut NamedTempFile) -> Result<()> {
        writeln!(file, "STAR Log.final.out")?;
        writeln!(file, "| Uniquely mapped reads % | 95.23% |")?;
        writeln!(file, "| Number of input reads | 1000000 |")?;
        writeln!(file, "| Uniquely mapped reads number | 952300 |")?;
        Ok(())
    }

    #[test]
    fn test_parse_star_log() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        write_valid_log(&mut temp_file)?;
        let stats = parse_star_log(temp_file.path().to_str().unwrap())?;
        assert_eq!(stats.mapping_ratio, "95.23%");
        assert_eq!(stats.total_reads, "1000000");
        assert_eq!(stats.uniquely_mapped_reads, "952300");
        assert!((stats.uniquely_mapped_ratio - 0.9523).abs() < 0.0001);
        Ok(())
    }

    #[test]
    fn test_parse_star_log_with_commas() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, "| Uniquely mapped reads % | 92.15% |")?;
        writeln!(temp_file, "| Number of input reads | 50,000,000 |")?;
        writeln!(temp_file, "| Uniquely mapped reads number | 46,075,000 |")?;
        let stats = parse_star_log(temp_file.path().to_str().unwrap())?;
        assert_eq!(stats.total_reads, "50000000");
        assert_eq!(stats.uniquely_mapped_reads, "46075000");
        Ok(())
    }

    #[test]
    fn test_missing_required_field_fails() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, "| Number of input reads | 1000000 |")?;
        writeln!(temp_file, "| Uniquely mapped reads number | 952300 |")?;
        assert!(parse_star_log(temp_file.path().to_str().unwrap()).is_err());
        Ok(())
    }

    #[test]
    fn duplicate_and_inconsistent_fields_are_rejected() -> Result<()> {
        let mut duplicate_file = NamedTempFile::new()?;
        write_valid_log(&mut duplicate_file)?;
        writeln!(duplicate_file, "| Number of input reads | 1000000 |")?;
        assert!(parse_star_log(duplicate_file.path().to_str().unwrap()).is_err());

        let mut inconsistent_file = NamedTempFile::new()?;
        writeln!(inconsistent_file, "| Uniquely mapped reads % | 50.0% |")?;
        writeln!(inconsistent_file, "| Number of input reads | 100 |")?;
        writeln!(inconsistent_file, "| Uniquely mapped reads number | 90 |")?;
        assert!(parse_star_log(inconsistent_file.path().to_str().unwrap()).is_err());
        Ok(())
    }
}
