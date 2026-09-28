# GV02-03 要求与交付文件

依据：[project-pool.pdf](course/project-pool.pdf) 第 20–22 页的 GV02-03，及其引用的代码仓库要求。本表核对技术内容与最终交付文件；课程成绩、实际答辩和真实过程记录由课程安排确定。

## 必做实验

| 原始要求 | 当前实现与证据 |
| --- | --- |
| 四维评价 pipeline | `src/gv_eval/`；1,000 条完整审计、136 条联合合格、114 条四维证据充分 |
| Pearson/Spearman、散点矩阵、冲突分析（RQ1） | 评价结果中的两个相关表；`reports/figures/metric_scatter_matrix.png`；项目报告第 4 节 |
| 同候选集上的至少三策略，比较重叠、质量和 Pareto 前沿大小（RQ2） | 三策略均在 114 条候选上生成 K=10/20/50；`strategy_summary.tsv`、`strategy_overlap.tsv`；`reports/pareto_front_summary.tsv` |
| 逐维消融与影响解释 | `ablation_summary.tsv`、`ablation_selections.json`；报告第 6 节 |
| 权重与阈值扰动、Jaccard/排名变化、敏感区域分析（RQ3） | 200 次 ±20% 权重扰动、27 组域/配对覆盖门槛，记录 Jaccard 和排序 Spearman；报告第 6 节 |
| 选做实验（最多两项） | 采用实验 7：1,000 次同池随机筛选；报告第 7 节。Pareto 图表用于解释必做策略比较，不另计下游功能验证 |

评价结果路径为 `results/gv02_03_v2/d_evaluation_member_a_v1/`。

## 最终交付

| 要求 | 文件 |
| --- | --- |
| 完整报告，回答 RQ1–RQ3 并给出适用场景 | [PROJECT_REPORT.md](../reports/PROJECT_REPORT.md) |
| 评价工具代码与筛选数据 | `src/`、`experiments/`、`results/` |
| 可复现仓库与 README 环境、运行说明 | [README](../README.md)、[复现说明](REPRODUCIBILITY.md)、`requirements-*.txt` |
| 答辩 PPT | [PROJECT_PRESENTATION.pptx](../reports/PROJECT_PRESENTATION.pptx) |
| 清晰的 data/preprocessing/models/experiments/results 结构 | 根目录同名文件夹 |

报告和 PPT 统一采用当前 1,000 条生成候选的批次。早期 M1–M5 过程文件与旧演示可在 `archive/pre-delivery-cleanup-20260928` 标签查看，不与当前报告混用。

## 完成范围

原文明确“本项目不训练新模型，而是构建评价和筛选系统”。因此，已发布候选的完整评价、比较与复现满足该任务的计算实验范围。GRU-VAE 训练是扩展，原 checkpoint 缺失限制原生成过程重放，但不阻断 GV02-03 的冻结候选评价。

原文没有要求通过湿实验确认功能。当前名单只能解释为计算筛选结果，参考仍为暂定版。家族建模和 Pfam/TM 扫描本轮核验发布证据，未将其描述为重新执行完成。

真实成员信息、会议、工时及实际答辩记录应来自真实活动；仓库整理不生成这些记录。原过程管理要求保存在 [process-requirements.docx](course/process-requirements.docx)。
