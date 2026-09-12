# 数据来源与使用说明

**候选入口：[candidates/t05/generated_gvp.fasta](candidates/t05/generated_gvp.fasta)，共 200 条。** 本项目围绕 T05 候选完成 GV02-03 四目标评价与筛选。计算结果提供候选优先级，但尚无湿实验核验的纯 GvpA 生成候选。

这里保留六份原始 FASTA，内容未经修改。`raw/` 是天然参考来源；文件名中的 GvpA 和头信息标签属于来源注释，不代表已验证家族或功能。

## 当前六份文件

| 文件 | 记录数 | 不同完整序列数 | 长度范围（aa） | 用途 |
| --- | ---: | ---: | ---: | --- |
| [candidates/t05/generated_gvp.fasta](candidates/t05/generated_gvp.fasta) | 200 | 200 | 66–512 | Transformer 生成候选，家族待核验 |
| [raw/t05/gvpa/GvpA_NotPartial.fasta](raw/t05/gvpa/GvpA_NotPartial.fasta) | 2000 | 1196 | 61–160 | T05 天然采集来源，保留注释与重复记录 |
| [raw/t05/gvpa/GvpA_RefSeq.fasta](raw/t05/gvpa/GvpA_RefSeq.fasta) | 862 | 862 | 54–177 | T05 RefSeq 来源 |
| [raw/t05/gvpa/rescued_GvpA_candidates.fasta](raw/t05/gvpa/rescued_GvpA_candidates.fasta) | 2 | 2 | 73–83 | 天然序列补救来源，名称中的 candidates 不表示模型生成 |
| [raw/design/GvpA.fasta](raw/design/GvpA.fasta) | 856 | 856 | 60–100 | 天然参考来源，在其余三份来源之外有 30 条独有序列 |
| [raw/t05/real_gvp.fasta](raw/t05/real_gvp.fasta) | 4379 | 4379 | 54–757 | 混合 Gvp 来源和原训练背景，用于来源核查及家族对照 |

“不同完整序列数”按序列字符串精确相同统计，不能替代按相似度聚类去冗余。四份名义 GvpA 来源的完整序列并集为 **1252 条**；按冻结规则清洗后正式参考集为 1224 条，再以 CD-HIT 90% identity 聚类得到 418 条 MSA/校准代表，均写入 `processed/gv02_03/`。

## 数据特点与评价边界

200 条候选均只含标准 20 种氨基酸，无空序列、重复 ID 或完全重复序列。长度均值 138.93 aa、中位数 103 aa，第一和第三四分位数为 93 和 114 aa。两条恰为生成上限 512 aa，应标记为可能截断，不能直接断定其截断或家族归属。

候选与四份名义 GvpA 来源并集、4379 条混合 Gvp 来源都没有完整序列完全相同的记录。这只排除了精确复制，尚未计算最大序列相似度、近重复簇或候选间多样性。当前无湿实验功能标签，不能据此计算功能预测准确率。

- `GvpA_NotPartial.fasta` 含大量重复序列记录，存在未知残基 `X`；NotPartial 命名不能替代质量核验。
- `raw/design/GvpA.fasta` 的头信息标注 GvpA 824 条、GvpJ 30 条、未标明上述家族 2 条，19 条含 `partial` 或 `fragment`。家族与完整性计数可重叠，不能将全部 856 条直接视为高质量 GvpA。
- 四份名义 GvpA 来源分别有其他三份未覆盖的完整序列：design 30、NotPartial 354、RefSeq 9、rescued 1。需保留来源与完整头信息，再进行可追溯的合并、去重和核验。
- design 独有的30条中，头信息标注 GvpA 3条、GvpJ 25条、未明确2条；保留价值包括潜在 GvpA 补充和注释溯源，不能把30条全部计作新增 GvpA。
- `real_gvp.fasta` 混合 GvpA、GvpC、GvpJ、GvpN 等家族，含 `X` 及不完整条目。其中头信息标注 GvpA 的 507 条序列已被四份来源覆盖；不能将全部 4379 条作为 GvpA 保守位点参考集。
- T05 候选只有生成编号，无亚家族标签。正式评分保留全部 200 条及无命中状态；37 条有 PF00741 可报告命中，27 条达到当前 GA+覆盖门槛。该结果只代表计算证据，不是家族或功能真值。

## 域标识修正

原任务记录中的 **PF01132 是 EFP（EF-P 的 OB 域），不是 GvpA 域**。相关的 **PF00741／Gas_vesicle** 家族也包含 GvpJ，因此命中只能支持相关家族证据，不能单独确认 GvpA 身份。本实验使用 PF00741.24，模型、来源元数据、域得分和覆盖度均已保存。[NCBI PF01132](https://www.ncbi.nlm.nih.gov/Structure/cdd/pfam01132)、[InterPro PF00741](https://www.ebi.ac.uk/interpro/entry/pfam/PF00741/)

## 来源与生成方法

T05 原代码和权重位于 [references/t05/src/model_training/](../references/t05/src/model_training/)，说明见 [model_training.md](../references/t05/doc/model_training.md)。生成设置为 temperature=0.8、top-k=10、最大长度 512、批量 200 条。缺失词表等依赖与 AutoDL 绝对路径尚未修复；当前未运行训练、补充生成或加载权重。

design 的 10 条旧候选和评分表已删除，不纳入候选池。design 天然 FASTA 因有独有序列而保留。原无效 MSA、结构预测输入和重复 FASTA 已删除；后续须从原始天然来源建立可靠比对。

最早来源路径映射见 [第一轮文件清单](../docs/cleanup/file_manifest.json)，第二轮删除记录见 [round2/deletions.json](../docs/cleanup/round2/deletions.json)。第一轮清单和 `fasta_inventory.json` 是历史记录，不能用作当前六份文件的清单。

## 更新输入审计与正式实验

在项目根目录执行：

```powershell
python -B -X utf8 preprocessing/analyze_inputs.py
```

该命令更新 [summary.json](../results/input_audit/summary.json) 与 [candidate_records.tsv](../results/input_audit/candidate_records.tsv)，包括文件数量、长度、SHA256、字符及注释审计和逐候选记录。它不清洗或覆盖原始 FASTA。

正式实验使用项目环境在仓库根目录运行：

```bash
python -B experiments/run_full_experiment.py --config configs/gv02_03.yaml
```

输出包括本目录下 `processed/gv02_03/` 的清洗参考、MSA、23 个统计保守位点和候选距离矩阵，以及 `results/gv02_03/` 的 HMMER/BLAST 原始证据、逐候选评分、策略名单、消融和鲁棒性结果。完整追溯信息见 `results/gv02_03/manifest.json`。
