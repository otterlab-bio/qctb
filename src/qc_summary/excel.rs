use crate::qc_summary::schema::{
    rnaseq_report_rows, standard_report_rows, ReportCell, ReportColumnSpec, ReportMode, ReportRow,
    REPORT_MISSING_VALUE, REPORT_SCHEMA_ID, REPORT_SCHEMA_NAME, REPORT_SCHEMA_VERSION,
};
use crate::qc_summary::{write_atomically, QCSummary, QCSummaryRNA};
use anyhow::{bail, Result};
use rust_xlsxwriter::{Color, Format, FormatAlign, Workbook, Worksheet};
use std::path::Path;

const MAX_EXACT_EXCEL_INTEGER: u64 = 9_007_199_254_740_992;

fn rounded(value: f64, decimal_places: Option<u8>) -> f64 {
    let factor = 10_f64.powi(i32::from(decimal_places.unwrap_or(0)));
    (value * factor).round() / factor
}

fn write_report_cell(
    sheet: &mut Worksheet,
    row: u32,
    column: u16,
    cell: &ReportCell,
    specification: &ReportColumnSpec,
    cell_format: &Format,
) -> Result<()> {
    match cell {
        ReportCell::Text(value) => {
            sheet.write_string_with_format(row, column, value, cell_format)?;
        }
        ReportCell::Integer(value) => {
            if *value > MAX_EXACT_EXCEL_INTEGER {
                bail!(
                    "Integer for column '{}' exceeds Excel's exact numeric range: {}",
                    specification.key,
                    value
                );
            }
            sheet.write_number_with_format(row, column, *value as f64, cell_format)?;
        }
        ReportCell::Decimal(value) => {
            sheet.write_number_with_format(
                row,
                column,
                rounded(*value, specification.decimal_places),
                cell_format,
            )?;
        }
        ReportCell::Missing => {
            sheet.write_string_with_format(row, column, REPORT_MISSING_VALUE, cell_format)?;
        }
    }
    Ok(())
}

fn add_report_sheet(workbook: &mut Workbook, rows: &[ReportRow], mode: ReportMode) -> Result<()> {
    let columns = mode.columns();
    let sheet = workbook.add_worksheet();
    sheet.set_name("Report")?;

    let header_format = Format::new()
        .set_bold()
        .set_font_color(Color::RGB(0xFFFFFF))
        .set_background_color(Color::RGB(0x4F81BD))
        .set_align(FormatAlign::Center);
    let cell_format = Format::new().set_align(FormatAlign::Left);

    for (column, specification) in columns.iter().enumerate() {
        sheet.write_string_with_format(
            0,
            column as u16,
            specification.excel_header,
            &header_format,
        )?;
        sheet.set_column_width(column as u16, 18)?;
    }

    for (row_index, row) in rows.iter().enumerate() {
        row.validate(mode)?;
        for (column_index, (cell, specification)) in row.cells.iter().zip(columns).enumerate() {
            write_report_cell(
                sheet,
                (row_index + 1) as u32,
                column_index as u16,
                cell,
                specification,
                &cell_format,
            )?;
        }
    }
    Ok(())
}

fn add_metadata_sheet(workbook: &mut Workbook, mode: ReportMode) -> Result<()> {
    let sheet = workbook.add_worksheet();
    sheet.set_name("qctb_metadata")?;
    sheet.write_string(0, 0, "schema_name")?;
    sheet.write_string(0, 1, REPORT_SCHEMA_NAME)?;
    sheet.write_string(1, 0, "schema_version")?;
    sheet.write_string(1, 1, REPORT_SCHEMA_VERSION)?;
    sheet.write_string(2, 0, "schema_id")?;
    sheet.write_string(2, 1, REPORT_SCHEMA_ID)?;
    sheet.write_string(3, 0, "mode")?;
    sheet.write_string(3, 1, mode.as_str())?;
    sheet.write_string(4, 0, "missing_value")?;
    sheet.write_string(4, 1, REPORT_MISSING_VALUE)?;

    let header_row = 6;
    for (column, header) in ["key", "excel_header", "type", "decimal_places"]
        .iter()
        .enumerate()
    {
        sheet.write_string(header_row, column as u16, *header)?;
    }
    for (row, specification) in mode.columns().iter().enumerate() {
        let row = header_row + 1 + row as u32;
        sheet.write_string(row, 0, specification.key)?;
        sheet.write_string(row, 1, specification.excel_header)?;
        sheet.write_string(row, 2, specification.column_type.as_str())?;
        if let Some(decimal_places) = specification.decimal_places {
            sheet.write_number(row, 3, f64::from(decimal_places))?;
        }
    }
    Ok(())
}

fn write_workbook(rows: &[ReportRow], output_path: &str, mode: ReportMode) -> Result<()> {
    let mut workbook = Workbook::new();
    add_report_sheet(&mut workbook, rows, mode)?;
    add_metadata_sheet(&mut workbook, mode)?;
    write_atomically(Path::new(output_path), |temporary_path| {
        workbook.save(temporary_path)?;
        Ok(())
    })
}

pub fn write_excel_standard_mode(
    summaries: &[QCSummary],
    output_path: &str,
    mode: ReportMode,
) -> Result<()> {
    let rows = standard_report_rows(summaries, mode)?;
    write_workbook(&rows, output_path, mode)
}

pub fn write_excel_standard(summaries: &[QCSummary], output_path: &str) -> Result<()> {
    write_excel_standard_mode(summaries, output_path, ReportMode::Standard)
}

pub fn write_excel_rnaseq(summaries: &[QCSummaryRNA], output_path: &str) -> Result<()> {
    let rows = rnaseq_report_rows(summaries)?;
    write_workbook(&rows, output_path, ReportMode::RnaSeq)
}
