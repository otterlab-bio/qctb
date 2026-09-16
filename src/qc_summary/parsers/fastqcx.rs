use super::stats::SeqkitStats;
use anyhow::{bail, Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

#[derive(Debug, Clone)]
pub struct FastqcxRow {
    pub num_seqs: u64,
    pub sum_len: u64,
    pub q20: f64,
    pub q30: f64,
    pub min_len: u32,
    pub avg_len: f64,
    pub max_len: u32,
}

fn required_value<'a>(
    values: &'a HashMap<String, String>,
    key: &str,
    file_path: &str,
) -> Result<&'a str> {
    values.get(key).map(String::as_str).with_context(|| {
        format!(
            "Key '{}' not found in Seqkit Statistics of {}",
            key, file_path
        )
    })
}

fn parse_u64_count(values: &HashMap<String, String>, key: &str, file_path: &str) -> Result<u64> {
    let raw_value = required_value(values, key, file_path)?;
    raw_value.parse::<u64>().with_context(|| {
        format!(
            "Field '{}' in {} must be an unsigned 64-bit integer, got '{}'",
            key, file_path, raw_value
        )
    })
}

fn parse_u32_count(values: &HashMap<String, String>, key: &str, file_path: &str) -> Result<u32> {
    let raw_value = required_value(values, key, file_path)?;
    raw_value.parse::<u32>().with_context(|| {
        format!(
            "Field '{}' in {} must be an unsigned 32-bit integer, got '{}'",
            key, file_path, raw_value
        )
    })
}

fn parse_finite_number(
    values: &HashMap<String, String>,
    key: &str,
    file_path: &str,
) -> Result<f64> {
    let raw_value = required_value(values, key, file_path)?;
    let parsed_value = raw_value.parse::<f64>().with_context(|| {
        format!(
            "Field '{}' in {} must be a number, got '{}'",
            key, file_path, raw_value
        )
    })?;
    if !parsed_value.is_finite() {
        bail!(
            "Field '{}' in {} must be finite, got '{}'",
            key,
            file_path,
            raw_value
        );
    }
    Ok(parsed_value)
}

fn parse_percentage(values: &HashMap<String, String>, key: &str, file_path: &str) -> Result<f64> {
    let percentage = parse_finite_number(values, key, file_path)?;
    if !(0.0..=100.0).contains(&percentage) {
        bail!(
            "Field '{}' in {} must be between 0 and 100, got '{}'",
            key,
            file_path,
            percentage
        );
    }
    Ok(percentage)
}

fn checked_add_counts(label: &str, left: u64, right: u64) -> Result<u64> {
    left.checked_add(right)
        .with_context(|| format!("Overflow while summing {}: {} + {}", label, left, right))
}

/// Parse the >>Seqkit Statistics section from a fastqc_data.txt file produced by fastqcx.
/// The section uses a tabular format: a header row (#file\tformat\t...) followed by
/// exactly one data row, with >>END_MODULE optionally appended to that row.
pub fn parse_fastqcx_data(file_path: &str) -> Result<FastqcxRow> {
    let file = File::open(file_path)
        .with_context(|| format!("Failed to open fastqcx data file: {}", file_path))?;
    let reader = BufReader::new(file);

    let mut in_section = false;
    let mut column_indices: HashMap<String, usize> = HashMap::new();
    let mut values_by_column: HashMap<String, String> = HashMap::new();
    let mut data_row_count = 0_u8;

    for line_result in reader.lines() {
        let line = line_result
            .with_context(|| format!("Failed to read fastqcx data file: {}", file_path))?;
        let trimmed_line = line.trim();
        if trimmed_line.starts_with(">>Seqkit Statistics") {
            in_section = true;
            continue;
        }
        if !in_section {
            continue;
        }

        let clean_line = trimmed_line.trim_end_matches(">>END_MODULE").trim_end();
        if clean_line.starts_with(">>") {
            break;
        }
        if clean_line.is_empty() {
            continue;
        }
        if clean_line.starts_with('#') {
            let headers = clean_line.trim_start_matches('#').split('\t');
            column_indices.clear();
            for (column_index, header) in headers.enumerate() {
                let normalized_header = header.trim().to_string();
                if column_indices
                    .insert(normalized_header.clone(), column_index)
                    .is_some()
                {
                    bail!(
                        "Duplicate Seqkit Statistics column '{}' in {}",
                        normalized_header,
                        file_path
                    );
                }
            }
            continue;
        }

        if column_indices.is_empty() {
            bail!(
                "Seqkit Statistics data row appears before its header in {}",
                file_path
            );
        }
        data_row_count = data_row_count.saturating_add(1);
        if data_row_count > 1 {
            bail!("Duplicate Seqkit Statistics data rows in {}", file_path);
        }

        let row_values: Vec<&str> = clean_line.split('\t').collect();
        for (column_name, column_index) in &column_indices {
            let raw_value = row_values.get(*column_index).with_context(|| {
                format!(
                    "Missing value for Seqkit Statistics column '{}' in {}",
                    column_name, file_path
                )
            })?;
            values_by_column.insert(column_name.clone(), raw_value.trim().to_string());
        }
    }

    if data_row_count == 0 {
        bail!("No Seqkit Statistics data row found in {}", file_path);
    }

    let average_length = parse_finite_number(&values_by_column, "avg_len", file_path)?;
    if average_length < 0.0 {
        bail!(
            "Field 'avg_len' in {} must be non-negative, got '{}'",
            file_path,
            average_length
        );
    }

    Ok(FastqcxRow {
        num_seqs: parse_u64_count(&values_by_column, "num_seqs", file_path)?,
        sum_len: parse_u64_count(&values_by_column, "sum_len", file_path)?,
        q20: parse_percentage(&values_by_column, "Q20(%)", file_path)?,
        q30: parse_percentage(&values_by_column, "Q30(%)", file_path)?,
        min_len: parse_u32_count(&values_by_column, "min_len", file_path)?,
        avg_len: average_length,
        max_len: parse_u32_count(&values_by_column, "max_len", file_path)?,
    })
}

/// Read 4 fastqcx fastqc_data.txt files and assemble into SeqkitStats
pub fn parse_seqkit_from_fastqcx(
    raw_r1: &str,
    raw_r2: &str,
    clean_r1: &str,
    clean_r2: &str,
) -> Result<SeqkitStats> {
    let r1 = parse_fastqcx_data(raw_r1)?;
    let r2 = parse_fastqcx_data(raw_r2)?;
    let c1 = parse_fastqcx_data(clean_r1)?;
    let c2 = parse_fastqcx_data(clean_r2)?;

    let reads_raw = checked_add_counts("raw reads", r1.num_seqs, r2.num_seqs)?;
    let bases_raw = checked_add_counts("raw bases", r1.sum_len, r2.sum_len)?;
    let reads_clean = checked_add_counts("clean reads", c1.num_seqs, c2.num_seqs)?;
    let bases_clean = checked_add_counts("clean bases", c1.sum_len, c2.sum_len)?;
    let clean_data_ratio = if bases_raw > 0 {
        bases_clean as f64 / bases_raw as f64
    } else {
        0.0
    };

    Ok(SeqkitStats {
        reads_raw_r1: r1.num_seqs,
        bases_raw_r1: r1.sum_len,
        q20_raw_r1: r1.q20,
        q30_raw_r1: r1.q30,
        min_len_raw_r1: r1.min_len,
        avg_len_raw_r1: r1.avg_len,
        max_len_raw_r1: r1.max_len,

        reads_raw_r2: r2.num_seqs,
        bases_raw_r2: r2.sum_len,
        q20_raw_r2: r2.q20,
        q30_raw_r2: r2.q30,
        min_len_raw_r2: r2.min_len,
        avg_len_raw_r2: r2.avg_len,
        max_len_raw_r2: r2.max_len,

        reads_clean_r1: c1.num_seqs,
        bases_clean_r1: c1.sum_len,
        q20_clean_r1: c1.q20,
        q30_clean_r1: c1.q30,
        min_len_clean_r1: c1.min_len,
        avg_len_clean_r1: c1.avg_len,
        max_len_clean_r1: c1.max_len,

        reads_clean_r2: c2.num_seqs,
        bases_clean_r2: c2.sum_len,
        q20_clean_r2: c2.q20,
        q30_clean_r2: c2.q30,
        min_len_clean_r2: c2.min_len,
        avg_len_clean_r2: c2.avg_len,
        max_len_clean_r2: c2.max_len,

        reads_raw,
        bases_raw,
        reads_clean,
        bases_clean,
        clean_data_ratio,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    struct TestFastqcxValues<'a> {
        num_seqs: &'a str,
        sum_len: &'a str,
        min_len: &'a str,
        avg_len: &'a str,
        max_len: &'a str,
        q20: &'a str,
        q30: &'a str,
        duplicate_data_row: bool,
    }

    impl Default for TestFastqcxValues<'static> {
        fn default() -> Self {
            Self {
                num_seqs: "100",
                sum_len: "15000",
                min_len: "50",
                avg_len: "150.0",
                max_len: "300",
                q20: "98.5",
                q30: "95.2",
                duplicate_data_row: false,
            }
        }
    }

    fn create_fastqcx_file(values: TestFastqcxValues<'_>) -> Result<NamedTempFile> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, ">>Seqkit Statistics\tpass")?;
        writeln!(
            temp_file,
            "#file\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tQ20(%)\tQ30(%)"
        )?;
        let data_row = format!(
            "test.fastq.gz\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            values.num_seqs,
            values.sum_len,
            values.min_len,
            values.avg_len,
            values.max_len,
            values.q20,
            values.q30
        );
        writeln!(temp_file, "{data_row}")?;
        if values.duplicate_data_row {
            writeln!(temp_file, "{data_row}")?;
        }
        writeln!(temp_file, ">>END_MODULE")?;
        Ok(temp_file)
    }

    #[test]
    fn test_parse_fastqcx_data() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, ">>Basic Statistics\tpass")?;
        writeln!(temp_file, "#Measure\tValue")?;
        writeln!(temp_file, "Filename\ttest.fastq.gz")?;
        writeln!(temp_file, ">>END_MODULE")?;
        writeln!(temp_file, ">>Seqkit Statistics\tpass")?;
        writeln!(temp_file, "#file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tQ1\tQ2\tQ3\tsum_gap\tN50\tN50_num\tQ20(%)\tQ30(%)\tAvgQual\tGC(%)\tsum_n")?;
        writeln!(temp_file, "test.fastq.gz\tFASTQ\tDNA\t1000000\t150000000\t50\t150.0\t300\t150\t150\t150\t0\t150\t1000000\t98.5\t95.2\t39.0\t50.0\t1000>>END_MODULE")?;

        let row = parse_fastqcx_data(temp_file.path().to_str().unwrap())?;
        assert_eq!(row.num_seqs, 1000000);
        assert_eq!(row.sum_len, 150000000);
        assert_eq!(row.q20, 98.5);
        assert_eq!(row.q30, 95.2);
        assert_eq!(row.min_len, 50);
        assert_eq!(row.avg_len, 150.0);
        assert_eq!(row.max_len, 300);

        Ok(())
    }

    #[test]
    fn count_fields_reject_invalid_unsigned_integers() -> Result<()> {
        for invalid_count in ["-1", "1.5", "NaN", "Inf", "18446744073709551616"] {
            let temp_file = create_fastqcx_file(TestFastqcxValues {
                num_seqs: invalid_count,
                ..TestFastqcxValues::default()
            })?;
            assert!(
                parse_fastqcx_data(temp_file.path().to_str().unwrap()).is_err(),
                "num_seqs value '{invalid_count}' should be rejected"
            );
        }

        let u32_overflow = create_fastqcx_file(TestFastqcxValues {
            min_len: "4294967296",
            ..TestFastqcxValues::default()
        })?;
        assert!(parse_fastqcx_data(u32_overflow.path().to_str().unwrap()).is_err());
        Ok(())
    }

    #[test]
    fn q20_and_q30_reject_non_finite_or_out_of_range_values() -> Result<()> {
        for invalid_percentage in ["-0.1", "100.1", "NaN", "Inf"] {
            let invalid_q20 = create_fastqcx_file(TestFastqcxValues {
                q20: invalid_percentage,
                ..TestFastqcxValues::default()
            })?;
            assert!(parse_fastqcx_data(invalid_q20.path().to_str().unwrap()).is_err());

            let invalid_q30 = create_fastqcx_file(TestFastqcxValues {
                q30: invalid_percentage,
                ..TestFastqcxValues::default()
            })?;
            assert!(parse_fastqcx_data(invalid_q30.path().to_str().unwrap()).is_err());
        }
        Ok(())
    }

    #[test]
    fn duplicate_seqkit_data_rows_are_rejected() -> Result<()> {
        let temp_file = create_fastqcx_file(TestFastqcxValues {
            duplicate_data_row: true,
            ..TestFastqcxValues::default()
        })?;
        let error = parse_fastqcx_data(temp_file.path().to_str().unwrap())
            .expect_err("duplicate data rows should fail");
        assert!(error
            .to_string()
            .contains("Duplicate Seqkit Statistics data rows"));
        Ok(())
    }

    #[test]
    fn summing_counts_rejects_overflow() {
        assert!(checked_add_counts("reads", u64::MAX, 1).is_err());
    }
}
