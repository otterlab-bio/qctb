use anyhow::{Context, Result};
use serde::Deserialize;
use serde_yaml::Value;
use std::collections::HashSet;
use std::path::{Path, PathBuf};

// ── New nested-format intermediate structs ─────────────────────────────────

#[derive(Debug, Deserialize, Default)]
struct NestedQCDirsQC {
    #[serde(default)]
    pub main: String,
    #[serde(default)]
    pub after: Option<String>,
    #[serde(default)]
    pub before: Option<String>,
}

#[derive(Debug, Deserialize, Default)]
struct NestedQCDirsBsmap {
    #[serde(default)]
    pub main: String,
}

#[derive(Debug, Deserialize, Default)]
struct NestedQCDirs {
    #[serde(default)]
    pub qc: NestedQCDirsQC,
    #[serde(default)]
    pub bsmap: NestedQCDirsBsmap,
    #[serde(default)]
    pub methylation_call: String,
    #[serde(default)]
    pub qualimap: String,
}

#[derive(Debug, Deserialize, Default)]
struct NestedOutput {
    #[serde(default)]
    pub trim_dir: String,
}

#[derive(Debug, Deserialize, Default)]
struct NestedWorkflowSpecies {
    #[serde(rename = "name", default)]
    _name: Vec<String>,
    #[serde(default)]
    graft: String,
}

#[derive(Debug, Deserialize, Default)]
struct NestedWorkflow {
    #[serde(default)]
    pub mode: String,
    #[serde(default)]
    pub species: NestedWorkflowSpecies,
}

#[derive(Debug, Deserialize, Default)]
struct NestedMetadata {
    #[serde(default)]
    pub sample_ids: Vec<String>,
}

// ── Raw config: handles both old flat and new nested format ────────────────

#[derive(Debug, Deserialize)]
struct RawConfig {
    // Old flat-format fields (all optional for backwards-compat)
    #[serde(rename = "SIDs", default)]
    sids: Vec<String>,
    #[serde(default)]
    graft: Option<String>,
    #[serde(default)]
    mode: Option<String>,
    #[serde(rename = "qcDir", default)]
    qc_dir: String,
    #[serde(rename = "trimDir", default)]
    trim_dir: String,
    #[serde(rename = "bsmapDir", default)]
    bsmap_dir: String,
    #[serde(rename = "outDir_mCall", default)]
    outdir_mcall: String,
    #[serde(rename = "outdir_qualimap", default)]
    qualimap_dir: String,
    #[serde(rename = "qcDir_before", default)]
    qcdir_before: Option<String>,
    #[serde(rename = "qcDir_after", default)]
    qcdir_after: Option<String>,

    // New nested-format fields
    #[serde(default)]
    metadata: NestedMetadata,
    #[serde(default)]
    directories: NestedQCDirs,
    #[serde(default)]
    output: NestedOutput,
    #[serde(default)]
    workflow: NestedWorkflow,
}

#[derive(Debug, Deserialize)]
struct CanonicalRunConfig {
    schema_version: String,
    run: CanonicalRunMetadata,
    workflow: CanonicalWorkflow,
    samples: Vec<CanonicalSample>,
    references: CanonicalReferences,
    paths: CanonicalPaths,
}

#[derive(Debug, Deserialize)]
struct CanonicalRunMetadata {
    immutable: bool,
}

#[derive(Debug, Deserialize)]
struct CanonicalWorkflow {
    scenario: String,
}

#[derive(Debug, Deserialize)]
struct CanonicalSample {
    id: String,
}

#[derive(Debug, Deserialize)]
struct CanonicalReferences {
    resolved: Vec<CanonicalReference>,
}

#[derive(Debug, Deserialize)]
struct CanonicalReference {
    role: String,
    id: String,
}

#[derive(Debug, Deserialize)]
struct CanonicalPaths {
    work: PathBuf,
}

// ── Public config returned to callers ─────────────────────────────────────

#[allow(non_snake_case)]
#[derive(Debug)]
pub struct QCConfig {
    pub SIDs: Vec<String>,
    pub graft: Option<String>,
    pub workflow_mode: Option<String>,
    pub qcDir: String,
    pub trimDir: String,
    pub bsmap_dir: String,
    pub qualimap_dir: String,
    pub outdir_mcall: String,
    pub qcdir_before: Option<String>,
    pub qcdir_after: Option<String>,
}

fn validate_sample_ids(sample_ids: &[String]) -> Result<()> {
    if sample_ids.is_empty() {
        anyhow::bail!("Configuration must contain at least one sample ID");
    }
    let mut seen = HashSet::new();
    for sample_id in sample_ids {
        if sample_id.trim().is_empty() {
            anyhow::bail!("Sample IDs must not be empty or whitespace-only");
        }
        if sample_id != sample_id.trim() {
            anyhow::bail!("Sample ID must not have leading or trailing whitespace: {sample_id:?}");
        }
        if sample_id.contains(['\t', '\r', '\n']) {
            anyhow::bail!("Sample ID contains a tab or newline: {sample_id:?}");
        }
        if !seen.insert(sample_id) {
            anyhow::bail!("Duplicate sample ID in configuration: '{sample_id}'");
        }
    }
    Ok(())
}

fn canonical_workflow_mode(scenario: &str) -> Result<&'static str> {
    match scenario {
        "rrbs" => Ok("RRBS"),
        "wgbs" => Ok("WGBS"),
        "bs-pdx" | "rna-pdx" => Ok("PDX"),
        "rnaseq" => Ok("RNASEQ"),
        _ => anyhow::bail!("Unsupported otter.run/v1 workflow scenario: {scenario}"),
    }
}

fn canonical_graft_reference(references: &[CanonicalReference]) -> Option<String> {
    references
        .iter()
        .find(|reference| reference.role == "graft")
        .or_else(|| {
            references
                .iter()
                .find(|reference| reference.role == "primary")
        })
        .map(|reference| reference.id.clone())
}

fn load_canonical_run_config(canonical: CanonicalRunConfig) -> Result<QCConfig> {
    if canonical.schema_version != "otter.run/v1" {
        anyhow::bail!(
            "Canonical QCTB configuration must use schema_version otter.run/v1, got {}",
            canonical.schema_version
        );
    }
    if !canonical.run.immutable {
        anyhow::bail!("Canonical QCTB configuration requires run.immutable: true");
    }
    if !canonical.paths.work.is_absolute() {
        anyhow::bail!("Canonical QCTB configuration requires an absolute paths.work");
    }

    let sample_ids = canonical
        .samples
        .into_iter()
        .map(|sample| sample.id)
        .collect::<Vec<_>>();
    validate_sample_ids(&sample_ids)?;
    let workflow_mode = canonical_workflow_mode(&canonical.workflow.scenario)?.to_owned();
    let graft = canonical_graft_reference(&canonical.references.resolved);
    if graft.is_none() {
        anyhow::bail!("Canonical QCTB configuration requires a graft or primary reference");
    }

    let work_directory = canonical.paths.work;
    let quality_control_directory = work_directory.join("QC");
    let methylation_directory = if workflow_mode == "RNASEQ" {
        work_directory.join("expression")
    } else {
        work_directory.join("mCall")
    };

    Ok(QCConfig {
        SIDs: sample_ids,
        graft,
        workflow_mode: Some(workflow_mode),
        qcDir: quality_control_directory.display().to_string(),
        trimDir: work_directory.join("trim").display().to_string(),
        bsmap_dir: work_directory.join("bsmap").display().to_string(),
        qualimap_dir: quality_control_directory
            .join("qualimap")
            .display()
            .to_string(),
        outdir_mcall: methylation_directory.display().to_string(),
        qcdir_before: Some(work_directory.join("fastqc_raw").display().to_string()),
        qcdir_after: Some(work_directory.join("fastqc_clean").display().to_string()),
    })
}

pub fn load_config(config_path: &Path) -> Result<QCConfig> {
    if config_path.is_dir() {
        anyhow::bail!(
            "QCTB --config requires a YAML file; use --config-dir for an Otter run directory: {}",
            config_path.display()
        );
    }
    let encoded = std::fs::read(config_path)
        .with_context(|| format!("Failed to open config file: {}", config_path.display()))?;
    let document: Value = serde_yaml::from_slice(&encoded)
        .with_context(|| format!("Failed to parse YAML config: {}", config_path.display()))?;

    if document
        .get("schema_version")
        .and_then(Value::as_str)
        .is_some_and(|schema_version| schema_version == "otter.run/v1")
    {
        let canonical: CanonicalRunConfig =
            serde_yaml::from_value(document).with_context(|| {
                format!(
                    "Failed to parse canonical run config: {}",
                    config_path.display()
                )
            })?;
        return load_canonical_run_config(canonical).with_context(|| {
            format!(
                "Failed to validate canonical run config: {}",
                config_path.display()
            )
        });
    }

    let raw: RawConfig = serde_yaml::from_value(document)
        .with_context(|| format!("Failed to parse YAML config: {}", config_path.display()))?;

    // Prefer old flat fields; fall back to new nested fields when flat fields are empty.
    let sids = if !raw.sids.is_empty() {
        raw.sids
    } else {
        raw.metadata.sample_ids
    };

    let qc_dir = if !raw.qc_dir.is_empty() {
        raw.qc_dir
    } else {
        raw.directories.qc.main.clone()
    };

    let trim_dir = if !raw.trim_dir.is_empty() {
        raw.trim_dir
    } else {
        raw.output.trim_dir.clone()
    };

    let bsmap_dir = if !raw.bsmap_dir.is_empty() {
        raw.bsmap_dir
    } else {
        raw.directories.bsmap.main.clone()
    };

    let qualimap_dir = if !raw.qualimap_dir.is_empty() {
        raw.qualimap_dir
    } else {
        raw.directories.qualimap.clone()
    };
    let outdir_mcall = if !raw.outdir_mcall.is_empty() {
        raw.outdir_mcall
    } else {
        raw.directories.methylation_call.clone()
    };

    let workflow_mode = raw.mode.or(if raw.workflow.mode.trim().is_empty() {
        None
    } else {
        Some(raw.workflow.mode)
    });
    let graft = raw.graft.or(if raw.workflow.species.graft.is_empty() {
        None
    } else {
        Some(raw.workflow.species.graft)
    });

    let qcdir_before = raw.qcdir_before.or(raw.directories.qc.before);
    let qcdir_after = raw.qcdir_after.or(raw.directories.qc.after);
    validate_sample_ids(&sids)?;

    Ok(QCConfig {
        SIDs: sids,
        graft,
        workflow_mode,
        qcDir: qc_dir,
        trimDir: trim_dir,
        bsmap_dir,
        qualimap_dir,
        outdir_mcall,
        qcdir_before,
        qcdir_after,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_load_config_old_format() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        writeln!(temp_file, "SIDs:")?;
        writeln!(temp_file, "  - sample1")?;
        writeln!(temp_file, "  - sample2")?;
        writeln!(temp_file, "graft: \"human\"")?;
        writeln!(temp_file, "qcDir: \"/qc\"")?;
        writeln!(temp_file, "trimDir: \"/trim\"")?;
        writeln!(temp_file, "bsmapDir: \"/bsmap\"")?;
        writeln!(temp_file, "outDir_mCall: \"/mcall\"")?;

        let config = load_config(temp_file.path())?;
        assert_eq!(config.SIDs.len(), 2);
        assert_eq!(config.SIDs[0], "sample1");
        assert_eq!(config.graft, Some("human".to_string()));
        assert_eq!(config.qcDir, "/qc");

        Ok(())
    }

    #[test]
    fn test_load_config_main_repository_nested_format() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        temp_file.write_all(
            br#"SIDs:
    - tumor_sample
    - normal_sample
directories:
    bsmap:
        main: /analysis/workflow/bsmap
    methylation_call: /analysis/workflow/mCall
    qualimap: /analysis/custom/qualimap
    qc:
        after: /analysis/workflow/fastqcx_clean
        before: /analysis/workflow/fastqcx_raw
        main: /analysis/workflow/QC
metadata:
    sample_ids:
        - tumor_sample
        - normal_sample
output:
    trim_dir: /analysis/workflow/trim
workflow:
    mode: PDX
    species:
        graft: human
        host: mouse
        name:
            - human
            - mouse
"#,
        )?;

        let config = load_config(temp_file.path())?;
        assert_eq!(
            config.SIDs,
            vec!["tumor_sample".to_string(), "normal_sample".to_string()]
        );
        assert_eq!(config.graft.as_deref(), Some("human"));
        assert_eq!(config.workflow_mode.as_deref(), Some("PDX"));
        assert_eq!(config.qcDir, "/analysis/workflow/QC");
        assert_eq!(
            config.qcdir_before.as_deref(),
            Some("/analysis/workflow/fastqcx_raw")
        );
        assert_eq!(
            config.qcdir_after.as_deref(),
            Some("/analysis/workflow/fastqcx_clean")
        );
        assert_eq!(config.trimDir, "/analysis/workflow/trim");
        assert_eq!(config.bsmap_dir, "/analysis/workflow/bsmap");
        assert_eq!(config.qualimap_dir, "/analysis/custom/qualimap");
        assert_eq!(config.outdir_mcall, "/analysis/workflow/mCall");

        Ok(())
    }

    #[test]
    fn test_load_config_from_immutable_otter_run_snapshot() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        temp_file.write_all(
            br#"schema_version: otter.run/v1
run:
  immutable: true
workflow:
  scenario: bs-pdx
samples:
  - id: tumor_sample
  - id: normal_sample
references:
  resolved:
    - role: graft
      id: hg38
    - role: host
      id: mm39
paths:
  work: /analysis/runs/run-example/work
"#,
        )?;

        let config = load_config(temp_file.path())?;
        assert_eq!(config.SIDs, vec!["tumor_sample", "normal_sample"]);
        assert_eq!(config.graft.as_deref(), Some("hg38"));
        assert_eq!(config.workflow_mode.as_deref(), Some("PDX"));
        assert_eq!(config.qcDir, "/analysis/runs/run-example/work/QC");
        assert_eq!(config.trimDir, "/analysis/runs/run-example/work/trim");
        assert_eq!(config.bsmap_dir, "/analysis/runs/run-example/work/bsmap");
        assert_eq!(
            config.qualimap_dir,
            "/analysis/runs/run-example/work/QC/qualimap"
        );
        assert_eq!(config.outdir_mcall, "/analysis/runs/run-example/work/mCall");

        Ok(())
    }

    #[test]
    fn immutable_otter_run_snapshot_is_required_for_canonical_loading() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        temp_file.write_all(
            br#"schema_version: otter.run/v1
run:
  immutable: false
workflow:
  scenario: rrbs
samples:
  - id: sample
references:
  resolved:
    - role: primary
      id: hg19
paths:
  work: /analysis/runs/run-example/work
"#,
        )?;

        let error = load_config(temp_file.path()).expect_err("mutable run must be rejected");
        assert!(
            error
                .chain()
                .any(|cause| cause.to_string().contains("run.immutable")),
            "unexpected canonical snapshot error: {error:#}"
        );
        Ok(())
    }

    #[test]
    fn empty_duplicate_and_whitespace_sample_ids_are_rejected() -> Result<()> {
        for yaml in [
            "SIDs: []\n",
            "SIDs:\n  - sample1\n  - sample1\n",
            "SIDs:\n  - ' sample1'\n",
            "SIDs:\n  - '   '\n",
        ] {
            let mut temp_file = NamedTempFile::new()?;
            temp_file.write_all(yaml.as_bytes())?;
            assert!(
                load_config(temp_file.path()).is_err(),
                "configuration should be rejected: {yaml:?}"
            );
        }
        Ok(())
    }
}
