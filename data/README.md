# 数据来源与使用说明

这里保留的是教师材料中的原始数据与候选数据，迁移未改动文件内容。文件名中的 GvpA 是来源命名，不代表已经完成家族和质量核验。

**当前候选入口：`candidates/t05/generated_gvp.fasta`。** 后续围绕 T05 的这200条生成序列识别 GvpA 候选并开展评价；`raw/` 中的天然序列用于参照，design 候选暂不并入正式候选池。目前还没有从T05中单独筛出并验证的纯GvpA生成候选文件。

## 文件清单

| 文件 | 记录数 | 不同序列数 | 长度范围（aa） | 用途 |
| --- | ---: | ---: | ---: | --- |
| `raw/design/GvpA.fasta` | 856 | 856 | 60–100 | 天然参考集来源之一 |
| `raw/t05/gvpa/GvpA_NotPartial.fasta` | 2000 | 1196 | 61–160 | T05 原始采集来源，保留头信息及重复条目 |
| `raw/t05/gvpa/GvpA_RefSeq.fasta` | 862 | 862 | 54–177 | T05 RefSeq 来源 |
| `raw/t05/gvpa/rescued_GvpA_candidates.fasta` | 2 | 2 | 73–83 | 天然序列补救来源，文件名中的 candidates 不代表模型生成 |
| `raw/t05/real_gvp.fasta` | 4379 | 4379 | 54–757 | 混合 Gvp 家族真实序列，保留已有训练基准及来源 |
| `candidates/t05/generated_gvp.fasta` | 200 | 200 | 66–512 | Transformer 生成的混合来源候选，尚未验证家族 |
| `candidates/design/candidates.fasta` | 10 | 10 | 83–99 | 原工程已经筛选过的小集合，可用于流程试跑 |

“不同序列数”按完整序列精确相同统计，不代表已经按序列相似度去冗余。

## 需要先处理的问题

- `raw/design/GvpA.fasta` 的头信息中，824 条标为 GvpA、30 条标为 GvpJ、2 条未明确标出上述家族；19 条含 `partial`。这些统计可能互相重叠，不能把 856 条直接视为高质量 GvpA。
- `GvpA_NotPartial.fasta` 有精确重复，部分序列含未知残基 `X`；名称里的 NotPartial 不能替代质量检查。`real_gvp.fasta` 也含 `X`。
- 四份名义 GvpA 来源按完整序列合并为 **1252 条不同序列**，这是尚未过滤的并集数量。四个文件分别都有其他三份未覆盖的信息，均予以保留，尚未生成正式合并参考集。
- T05 的 `real_gvp.fasta` 包含多种 Gvp 家族。其中头信息标有 GvpA 的 507 条已被上述四份来源覆盖；整个混合集有 532 条序列与四份来源并集完全一致。头信息只提供注释线索，不是家族确认结果。
- T05 的 200 条候选没有家族标签。应先明确 GvpA 候选池的构建和约束评价方案，不能把全体直接当作已通过家族核验的 GvpA。
- design 的 10 条候选由 8 条 Transformer、2 条 Diffusion 序列组成，是已筛选子集。完整生成池和模型权重不在材料中；只在这 10 条上比较策略会受到此前筛选的影响。
- design 候选与 T05 候选之间没有完全相同的序列；这不等同于已证明它们具有高新颖性或高多样性。

## 来源与副本处理

原始路径到新路径的完整映射见 [file_manifest.json](../docs/cleanup/file_manifest.json)。

- `design/data/natural_gvp.fasta` 与保留的 `raw/design/GvpA.fasta` 字节完全一致，已删除副本。
- T05 的 `deduplicated.fasta`、A-Prot 输入中的 `real_gvp.fasta` 与保留的 `raw/t05/real_gvp.fasta` 字节完全一致，已删除副本。
- A-Prot 输入中的 `generated_gvp.fasta` 与保留的 T05 候选字节完全一致，已删除副本。
- 原 `design/data/*_msa.a3m` 全部由原序列截成统一长度产生。GvpA 文件把 856 条序列全部截至 60 aa，不能作为可靠的 MSA，已经删除。后续应从保留的原始序列建立真正的比对。
- T05 的混合家族结构预测输入及其衍生对齐已删除。它们不应作为 GvpA 关键位点的参考比对。

## 生成方法与旧分数

T05 的来源代码与已有权重位于 `../references/t05/src/model_training/`，原文说明位于 `../references/t05/doc/model_training.md`。脚本使用 temperature=0.8、top-k=10、最大长度512、生成200条；缺失依赖和绝对路径尚未修复，未验证重新生成能否复现现存文件。

design 的候选来源表位于 `../references/design/candidate_metadata/candidates_detailed.csv`。其中旧分数使用另一套五维体系，包含启发式计算，保留用于追溯，不作为本项目四维指标。详细限制见 [参考说明](../references/README.md)。

可从根目录运行 `python -X utf8 preprocessing/audit_fasta.py` 重新核对单文件统计。完整快照见 [fasta_inventory.json](../docs/cleanup/fasta_inventory.json)。
