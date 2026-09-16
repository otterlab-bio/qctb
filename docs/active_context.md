# System Context (Updated: 2026-07-25)

## 1. 已实现的核心模块

### Configuration Loading
- **Path**: `src/qc_summary/config.rs`
- **Contract**: 同时加载旧平铺配置和 otter 嵌套 YAML；保留 `workflow.mode`、`workflow.species.graft` 和 `directories.qualimap`。
- **Validation**: SID 列表不能为空；拒绝空白 SID、重复 SID、首尾空白及 tab/newline。

### QC Aggregation
- **Path**: `src/qc_summary/aggregator.rs`
- **Public Methods**: `process_all_samples`, `process_all_samples_rnaseq`。
- **fastqcx Contract**: 从 `<sample>_*_fastqcx/fastqc_data.txt` 读取 FastQC-compatible Seqkit Statistics，不提供旧产品路径或符号别名。
- **methx Contract**: 从 `methylation_call/methrixh5` 读取 Methrix QC 产物；报告缺失为可选，报告存在但目标样本缺失或重复时失败。

### Versioned Native Report
- **Path**: `src/qc_summary/schema.rs`, `src/qc_summary/tsv.rs`, `src/qc_summary/excel.rs`
- **Schema**: `qctb.report/1.0.0`。
- **Modes**: RRBS、WGBS、RNA-seq、PDX；PDX v1 与标准 BS 列契约相同，模式身份单独记录。
- **Contract**: 统一固定 Excel/TSV 列名、顺序、类型、十进制位数和 `N/A` 缺失值。
- **Metadata**: TSV 前两行为 schema/mode 注释；Excel 使用 `qctb_metadata` 工作表记录 schema、mode 与逐列定义。
- **Publication**: Excel 和 TSV 均经同目录临时文件、`sync_all` 和原子替换发布。

### Strict Parsers
- **fastqcx**: count 使用精确整数解析；拒绝负数、小数、NaN、Inf、溢出、重复数据行和越界 Q20/Q30。
- **Bismark/STAR**: 必需字段必须且只能出现一次；百分比范围、非零分母、映射数不超过总数，并校验报告百分比与计数一致。
- **Qualimap**: 重复字段失败；mapping quality、duplicated reads 与 duplication rate 严格解析并保留 `%` 单位。
- **Trim Galore**: 字段按完整行匹配，重复或缺失字段失败。
- **methx**: Methrix workbook 在测试中运行时生成；校验唯一列名、必需指标、唯一样本、整数 count、percent 范围及 coverage 阈值单调性。

## 2. 关键类型

| Type | Path | Purpose |
|---|---|---|
| `QCConfig` | `src/qc_summary/config.rs` | `SIDs`, `workflow_mode`, `graft`, fastqcx/trim/bsmap/Qualimap/Methrix paths |
| `QCSummary` / `QCSummaryRNA` | `src/qc_summary/aggregator.rs` | BS/RNA 样本汇总 |
| `FastqcxRow` | `src/qc_summary/parsers/fastqcx.rs` | fastqcx 的单个 FastQC-compatible Seqkit Statistics 数据行 |
| `SeqkitStats` | `src/qc_summary/parsers/stats.rs` | fastqcx Seqkit Statistics 聚合类型 |
| `MethrixCoverageRow` / `MethrixAnnotationBySampleRow` | `src/qc_summary/parsers/methx.rs` | methx 输出的 Methrix 科学数据行 |
| `ReportMode` | `src/qc_summary/schema.rs` | RRBS/WGBS/RNA-seq/PDX schema mode |
| `ReportColumnSpec` | `src/qc_summary/schema.rs` | 列 key、Excel header、类型和舍入 |
| `ReportCell` / `ReportRow` | `src/qc_summary/schema.rs` | writer 共享 typed output |

## 3. Golden Coverage

- `tests/golden/rrbs.tsv.golden`
- `tests/golden/wgbs.tsv.golden`
- `tests/golden/rnaseq.tsv.golden`
- `tests/golden/pdx.tsv.golden`
- `tests/report_schema.rs` 同时验证 Excel schema metadata 与列数。

## 4. 验证状态

- `Cargo.lock` 已生成并由应用仓跟踪。
- Rust fmt/check/clippy/test 门禁使用 `--locked`。
- 当前本地验证：`cargo fmt --check`、`cargo test --locked`（40 library、1 CLI、6 integration）、`cargo clippy --locked --all-targets --all-features -- -D warnings` 全部通过。

## 5. 剩余外部集成

- [ ] 使用 otter 真实 Snakemake 配置分别运行 RRBS/WGBS/RNA/PDX 数据链。
- [ ] 与 `methx` 实际生成的 coverage/annotation workbook 做跨仓 smoke test。
- [ ] PDX v2 若新增 host/graft 双物种指标，必须升级 schema 版本，不得静默扩展 v1 列。
