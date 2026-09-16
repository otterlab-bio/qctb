pub mod aggregator;
pub mod config;
pub mod excel;
pub mod parsers;
pub mod schema;
pub mod tsv;

use anyhow::{Context, Result};
use std::path::Path;
use tempfile::Builder;

pub fn write_atomically<F>(output_path: &Path, write_temporary_output: F) -> Result<()>
where
    F: FnOnce(&Path) -> Result<()>,
{
    let parent_directory = output_path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    std::fs::create_dir_all(parent_directory).with_context(|| {
        format!(
            "Failed to create output directory: {}",
            parent_directory.display()
        )
    })?;

    let temporary_suffix = output_path
        .extension()
        .and_then(|extension| extension.to_str())
        .map(|extension| format!(".{extension}"))
        .unwrap_or_else(|| ".tmp".to_string());
    let temporary_file = Builder::new()
        .prefix(".qctb-")
        .suffix(&temporary_suffix)
        .tempfile_in(parent_directory)
        .with_context(|| {
            format!(
                "Failed to create temporary output in {}",
                parent_directory.display()
            )
        })?;

    write_temporary_output(temporary_file.path())?;
    std::fs::File::open(temporary_file.path())
        .with_context(|| {
            format!(
                "Failed to reopen temporary output: {}",
                temporary_file.path().display()
            )
        })?
        .sync_all()
        .with_context(|| {
            format!(
                "Failed to synchronize temporary output: {}",
                temporary_file.path().display()
            )
        })?;

    temporary_file
        .persist(output_path)
        .map_err(|persist_error| {
            anyhow::anyhow!(
                "Failed to atomically replace output {}: {}",
                output_path.display(),
                persist_error.error
            )
        })?;
    Ok(())
}

pub use aggregator::{
    process_all_samples, process_all_samples_rnaseq, process_sample, QCSummary, QCSummaryRNA,
};
pub use config::{load_config, QCConfig};
pub use excel::{write_excel_rnaseq, write_excel_standard, write_excel_standard_mode};
pub use schema::{ReportMode, REPORT_SCHEMA_ID, REPORT_SCHEMA_NAME, REPORT_SCHEMA_VERSION};
pub use tsv::{write_tsv_rnaseq, write_tsv_standard};

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::bail;
    use tempfile::TempDir;

    #[test]
    fn atomic_write_creates_parent_directories() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let output_path = temp_dir.path().join("nested").join("summary.tsv");

        write_atomically(&output_path, |temporary_path| {
            std::fs::write(temporary_path, "complete output")?;
            Ok(())
        })?;

        assert_eq!(std::fs::read_to_string(output_path)?, "complete output");
        Ok(())
    }

    #[test]
    fn atomic_write_failure_preserves_existing_output() -> Result<()> {
        let temp_dir = TempDir::new()?;
        let output_path = temp_dir.path().join("summary.tsv");
        std::fs::write(&output_path, "previous complete output")?;

        let result = write_atomically(&output_path, |temporary_path| {
            std::fs::write(temporary_path, "partial replacement")?;
            bail!("simulated output failure")
        });

        assert!(result.is_err());
        assert_eq!(
            std::fs::read_to_string(output_path)?,
            "previous complete output"
        );
        Ok(())
    }
}
