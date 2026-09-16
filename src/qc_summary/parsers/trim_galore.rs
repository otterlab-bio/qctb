use anyhow::{Context, Result};
use regex::Regex;
use std::fs::read_to_string;

#[derive(Debug, Clone)]
pub struct TrimStats {
    // R1 stats
    pub reads_with_adapter_r1: String,
    pub reads_write_r1: String,
    pub bp_qc_remove_r1: String,
    pub bp_write_r1: String,

    // R2 stats
    pub reads_with_adapter_r2: String,
    pub reads_write_r2: String,
    pub bp_qc_remove_r2: String,
    pub bp_write_r2: String,
}

pub fn parse_trim_report(file_path: &str) -> Result<TrimStats> {
    let content = read_to_string(file_path)
        .with_context(|| format!("Failed to read trim report: {}", file_path))?;

    let re_adapter = Regex::new(r"(?m)^\s*Reads with adapters:\s*(\S.*)\s*$")?;
    let re_written = Regex::new(r"(?m)^\s*Reads written \(passing filters\):\s*(\S.*)\s*$")?;
    let re_qc = Regex::new(r"(?m)^\s*Quality-trimmed:\s*(\S.*)\s*$")?;
    let re_total = Regex::new(r"(?m)^\s*Total written \(filtered\):\s*(\S.*)\s*$")?;

    let extract = |regex: &Regex, field: &str| -> Result<String> {
        let matches = regex
            .captures_iter(&content)
            .filter_map(|captures| captures.get(1))
            .map(|value| value.as_str().trim().to_string())
            .collect::<Vec<_>>();
        match matches.as_slice() {
            [value] => Ok(value.clone()),
            [] => anyhow::bail!(
                "Missing required Trim Galore field '{}' in {}",
                field,
                file_path
            ),
            _ => anyhow::bail!(
                "Duplicate Trim Galore field '{}' appears {} times in {}",
                field,
                matches.len(),
                file_path
            ),
        }
    };

    Ok(TrimStats {
        reads_with_adapter_r1: extract(&re_adapter, "Reads with adapters")?,
        reads_write_r1: extract(&re_written, "Reads written (passing filters)")?,
        bp_qc_remove_r1: extract(&re_qc, "Quality-trimmed")?,
        bp_write_r1: extract(&re_total, "Total written (filtered)")?,
        // Note: R2 stats will be filled when parsing R2 file
        reads_with_adapter_r2: String::new(),
        reads_write_r2: String::new(),
        bp_qc_remove_r2: String::new(),
        bp_write_r2: String::new(),
    })
}

pub fn parse_trim_reports(r1_file: &str, r2_file: &str) -> Result<TrimStats> {
    let r1_stats = parse_trim_report(r1_file)?;
    let r2_stats = parse_trim_report(r2_file)?;

    Ok(TrimStats {
        reads_with_adapter_r1: r1_stats.reads_with_adapter_r1,
        reads_write_r1: r1_stats.reads_write_r1,
        bp_qc_remove_r1: r1_stats.bp_qc_remove_r1,
        bp_write_r1: r1_stats.bp_write_r1,
        // parse_trim_report always populates _r1 fields; use those for the R2 file's data
        reads_with_adapter_r2: r2_stats.reads_with_adapter_r1,
        reads_write_r2: r2_stats.reads_write_r1,
        bp_qc_remove_r2: r2_stats.bp_qc_remove_r1,
        bp_write_r2: r2_stats.bp_write_r1,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_parse_trim_report() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, "This is a trim galore report")?;
        writeln!(temp_file, "Reads with adapters:        50000 (5.0%)")?;
        writeln!(temp_file, "Reads written (passing filters): 950000")?;
        writeln!(temp_file, "Quality-trimmed: 2500000 bp")?;
        writeln!(temp_file, "Total written (filtered): 142500000 bp")?;

        let stats = parse_trim_report(temp_file.path().to_str().unwrap())?;
        assert_eq!(stats.reads_with_adapter_r1, "50000 (5.0%)");
        assert_eq!(stats.reads_write_r1, "950000");

        Ok(())
    }

    #[test]
    fn test_parse_trim_reports() -> Result<()> {
        // Create R1 file
        let mut r1_file = NamedTempFile::new()?;
        writeln!(r1_file, "Reads with adapters:        50000")?;
        writeln!(r1_file, "Reads written (passing filters): 950000")?;
        writeln!(r1_file, "Quality-trimmed: 2500000 bp")?;
        writeln!(r1_file, "Total written (filtered): 142500000 bp")?;

        // Create R2 file
        let mut r2_file = NamedTempFile::new()?;
        writeln!(r2_file, "Reads with adapters:        60000")?;
        writeln!(r2_file, "Reads written (passing filters): 940000")?;
        writeln!(r2_file, "Quality-trimmed: 2600000 bp")?;
        writeln!(r2_file, "Total written (filtered): 141000000 bp")?;

        let stats = parse_trim_reports(
            r1_file.path().to_str().unwrap(),
            r2_file.path().to_str().unwrap(),
        )?;

        assert_eq!(stats.reads_with_adapter_r1, "50000");
        assert_eq!(stats.reads_with_adapter_r2, "60000");

        Ok(())
    }

    #[test]
    fn duplicate_fields_fail() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, "Reads with adapters: 50000")?;
        writeln!(temp_file, "Reads with adapters: 60000")?;
        writeln!(temp_file, "Reads written (passing filters): 950000")?;
        writeln!(temp_file, "Quality-trimmed: 2500000 bp")?;
        writeln!(temp_file, "Total written (filtered): 142500000 bp")?;
        assert!(parse_trim_report(temp_file.path().to_str().unwrap()).is_err());
        Ok(())
    }

    #[test]
    fn test_missing_required_field_fails() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, "Reads with adapters: 50000")?;
        assert!(parse_trim_report(temp_file.path().to_str().unwrap()).is_err());
        Ok(())
    }
}
