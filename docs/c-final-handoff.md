# C→下游：真实 VAE 批次统一复核交接（member A v1）

本交接把 C 已发布的主 QC/相似度、Pfam 域架构和跨膜预测结果按候选 ID 汇总，供后续分析和人工复核使用。它是**预警性证据索引**，不是 GvpA 身份、结构或功能的最终验证，也不修改 B 的候选池、基础 QC 或家族标签。

## 输入版本与可追溯性

- B 冻结批次：`origin/feat/sequence-vae` 的 `10d9da8`（`member_a_v1_seed42`），1,000 条候选、346 条暂定参考；主分析池 147、歧义探索池 47、排除池 806。
- C 主检查：[原始交接报告](../results/gv02_03_c_handoff_member_a_v1/report.md)与 `candidate_audit.tsv`，覆盖导出元数据/EOS、基础 QC、序列模式告警、全局最近参考匹配及独立距离矩阵。
- 后续独立检查：[Pfam 域架构](../results/gv02_03_c_domains_member_a_v1/README.md)及[跨膜拓扑](../results/gv02_03_c_tm_member_a_v1/README.md)。原始报告中“尚未检查”只描述**当时那次运行**；这两项是在同一冻结输入上随后补充的结果，未回写历史产物。

汇总脚本先逐一核对三批已发布结果清单的输出 SHA-256，并比较候选、参考和分池输入哈希；对 1,000 条候选按 ID 检查序列 SHA-256、长度和原池标签一致，拒绝漏行、重复 ID 或字段冲突。它还依据各自经哈希核验的参考集阈值检查 Pfam 命中数、警告状态和跨膜片段数，来源表内部矛盾时不生成交接结果。新结果的输入/输出哈希及汇总脚本规范化换行后的 SHA-256 见[汇总清单](../results/gv02_03_c_final_handoff_member_a_v1/manifest.json)。

## 交付文件与判读

- [`candidate_c_review.tsv`](../results/gv02_03_c_final_handoff_member_a_v1/candidate_c_review.tsv)：全部 1,000 条，保留原 `candidate_audit.tsv` 字段和顺序，追加域及跨膜状态、坐标/片段数、`c_review_reasons` 与 `c_manual_review`。
- [`manual_review_c.tsv`](../results/gv02_03_c_final_handoff_member_a_v1/manual_review_c.tsv)：上述总表中 `c_manual_review=True` 的原顺序子集，915 条。
- [`summary.json`](../results/gv02_03_c_final_handoff_member_a_v1/summary.json) 与 `manifest.json`：计数、边界和哈希。

原 `manual_review` 列保持 C 主检查的判定；新 `c_manual_review` 是“原复核名单 ∪ Pfam 多段/非目标警告 ∪ 无显著 Pfam 命中的未解析状态 ∪ 跨膜预测警告”。`c_review_reasons` 用 `base:`、`domain:`、`tm:` 前缀标明来源。它是**人工复核标记，不是淘汰或通过证明**。尤其 `domain:no_pfam_hit_unresolved` 表示未取得显著 Pfam 命中，不能解释为“无结构域”。

## 本批结果

| 原候选池 | 全部候选 | 合并后人工复核 |
| --- | ---: | ---: |
| `main_supported_gvpa` | 147 | 62 |
| `ambiguous_exploration` | 47 | 47 |
| `excluded` | 806 | 806 |
| **合计** | **1,000** | **915** |

Pfam 扫描有 **30 条多段 PF00741 命中**需要复核；跨膜预测有 **7 条警告**，两组不重叠。另有 **395 条无显著 Pfam 命中**，均在原排除池，属于证据未解析。37 条新增的结构域/跨膜警告 ID 早已在原 915 条复核名单中，所以合并后仍为 915，**不能将 915 与 37 相加**。本批未检出非 PF00741 的显著 Pfam 命中；没有检测到不等于证明不存在。逐段证据请回看两项独立检查的命中表/预测表。

## 给下游的使用边界

1. 主分析池与歧义池、排除池继续分开；`qc_pass`、`c_manual_review=False` 或 Pfam 命中都不自动赋予 GvpA 身份。不要从排除池补足主名单。
2. 最近参考是全局比对证据，不代替 B 的竞争性家族鉴定；未解析距离保留缺失，不能补成最大距离。参考集仍为暂定版本。
3. 多段 PF00741 可能是拆分或重叠比对，不能直接称为多份完整结构域；跨膜标签是计算预测，不是实验证据。对 30 条和 7 条警告应查看逐条坐标并人工判断，暂不据此重分池。
4. EOS/生成上限核验来自导出元数据和序列一致性，未重放 B checkpoint 或原始 token 流。这里没有最终候选推荐、独立结构预测或功能/表达/组装验证。

## 本地复现

在仓库根目录、三个已发布 C 结果目录均存在时运行（使用空的新输出目录）：

```powershell
$env:PYTHONPATH = (Resolve-Path src).Path
.venv/Scripts/python.exe -m gv_eval.c_final_handoff `
  --base results/gv02_03_c_handoff_member_a_v1 `
  --domains results/gv02_03_c_domains_member_a_v1 `
  --tm results/gv02_03_c_tm_member_a_v1 `
  --output results/gv02_03_c_final_handoff_member_a_v1_rerun
```

这一步只读取已发布的结果文件，不需要重新下载 Pfam 或运行 PureseqTM；脚本不会覆盖已有输出。三个原结果如何生成及其工具限制，见各自的说明与 manifest。
