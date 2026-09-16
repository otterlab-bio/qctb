use crate::qc_summary::schema::{
    rnaseq_report_rows, standard_report_rows, ReportMode, ReportRow, REPORT_SCHEMA_ID,
};
use crate::qc_summary::{write_atomically, QCSummary, QCSummaryRNA};
use anyhow::Result;
use std::io::Write;
use std::path::Path;

fn write_tsv_contents(rows: &[ReportRow], mode: ReportMode, output_path: &Path) -> Result<()> {
    let columns = mode.columns();
    let mut file = std::fs::File::create(output_path)?;
    writeln!(file, "# qctb_schema={REPORT_SCHEMA_ID}")?;
    writeln!(file, "# qctb_mode={}", mode.as_str())?;
    writeln!(
        file,
        "{}",
        columns
            .iter()
            .map(|column| column.key)
            .collect::<Vec<_>>()
            .join("\t")
    )?;

    for row in rows {
        row.validate(mode)?;
        writeln!(
            file,
            "{}",
            row.cells
                .iter()
                .zip(columns)
                .map(|(cell, column)| cell.tsv_value(column))
                .collect::<Result<Vec<_>>>()?
                .join("\t")
        )?;
    }
    Ok(())
}

pub fn write_tsv_standard(
    summaries: &[QCSummary],
    output_path: &Path,
    mode: ReportMode,
) -> Result<()> {
    let rows = standard_report_rows(summaries, mode)?;
    write_atomically(output_path, |temporary_path| {
        write_tsv_contents(&rows, mode, temporary_path)
    })
}

pub fn write_tsv_rnaseq(summaries: &[QCSummaryRNA], output_path: &Path) -> Result<()> {
    let rows = rnaseq_report_rows(summaries)?;
    write_atomically(output_path, |temporary_path| {
        write_tsv_contents(&rows, ReportMode::RnaSeq, temporary_path)
    })
}
