# GV02-03 要求与交付文件

依据：[project-pool.pdf](course/project-pool.pdf) 第 20–22 页的 GV02-03，及其引用的代码仓库要求。2026-10-07 核对结果：四项必做及选做实验 6、7 已完成，选做实验 5 未做。本表核对技术内容与交付文件；课程成绩、实际答辩和真实过程记录由课程安排确定。

## 必做实验

| 原始要求 | 当前实现与证据 |
| --- | --- |
| 四维评价 pipeline | `src/gv_eval/`；1,000 条完整审计、136 条联合合格、114 条四维证据充分 |
| Pearson/Spearman、散点矩阵、冲突分析（RQ1） | 评价结果中的两个相关表；`reports/figures/metric_scatter_matrix.png`；项目报告第 4 节 |
| 同候选集上的至少三策略，比较重叠、质量和 Pareto 前沿大小（RQ2） | 三策略均在 114 条候选上生成 K=10/20/50；`strategy_summary.tsv`、`strategy_overlap.tsv`；`reports/pareto_front_summary.tsv` |
| 逐维消融与影响解释 | `ablation_summary.tsv`、`ablation_selections.json`；报告第 6 节 |
| 权重与阈值扰动、Jaccard/排名变化、敏感区域分析（RQ3） | 200 次 ±20% 权重扰动、27 组域/配对覆盖门槛，记录 Jaccard 和排序 Spearman；报告第 6 节 |

评价结果路径为 `results/gv02_03_v2/d_evaluation_member_a_v1/`。

## 选做实验：完成两项

原文第 21 页要求“从中至多选 2 项”。当前选择如下：

| 实验 | 状态 | 实现与证据 |
| --- | --- | --- |
| 5：下游验证代理 | 未做 | 未接入 GV02-02 性质预测模型；现有域和跨膜检查不计作此项 |
| 6：Pareto 前沿分析 | 已完成 | 四维平行坐标及逐候选优势；10 条前沿与 104 条非前沿的六项预定特征、完整 20 种氨基酸组成；一对一长度匹配及 100 次等价最优匹配敏感性检查。见[实验报告](../results/gv02_03_v2/pareto_analysis_member_a_v1/EXPERIMENT_6_REPORT.md) |
| 7：与随机筛选对比 | 已完成 | 同一 114 条池、相同 Top-20 预算，1,000 次随机抽样，比较各目标及距离覆盖；见项目报告第 7 节 |

实验 6 的 100 次匹配保留同一个最小总长度差，只检查等价对照选择的影响；它们不是独立生物学重复。其特征分析不改变原始四维目标与候选资格，也不另计为新的选做项目。

## 最终交付

| 要求 | 文件 |
| --- | --- |
| 完整报告，回答 RQ1–RQ3 并给出适用场景 | [PROJECT_REPORT.md](../reports/PROJECT_REPORT.md)及[实验 6 专项报告](../results/gv02_03_v2/pareto_analysis_member_a_v1/EXPERIMENT_6_REPORT.md) |
| 评价工具代码与筛选数据 | `src/`、`experiments/`、`results/` |
| 可复现仓库与 README 环境、运行说明 | [README](../README.md)、[复现说明](REPRODUCIBILITY.md)、`requirements-*.txt` |
| 答辩 PPT | [PROJECT_PRESENTATION.pptx](../reports/PROJECT_PRESENTATION.pptx) |
| 清晰的 data/preprocessing/models/experiments/results 结构 | 根目录同名文件夹 |

报告和 PPT 均基于当前 1,000 条生成候选的批次；新补充的实验 6 已进入报告，已有演示稿尚未同步该扩展。早期 M1–M5 过程文件与旧演示可在 `archive/pre-delivery-cleanup-20260928` 标签查看，不与当前报告混用。

## 评分与过程材料

[课程介绍](course/course-introduction.pdf)第 9 页标明总评由 **40% 平时过程、60% 项目成果**构成。观察点包括问题建模、实验设计与分析、方法选择依据、可复现性、真实团队协作及答辩表达；成果部分子项标注与总项的折算口径存在歧义，细分口径以教师为准，不据此自行估算得分。

课程介绍第 18–20 页要求研究问题驱动、真实 Git 与在线协作记录、关键实验命令行复现。补齐实验可以加强技术交付，不能替代真实会议、参与和答辩。课程介绍按原始字节从 `archive/pre-delivery-cleanup-20260928:intro-mlproj.pdf` 恢复保存；会议与过程要求另见 [process-requirements.docx](course/process-requirements.docx)。

## 完成范围

原文明确“本项目不训练新模型，而是构建评价和筛选系统”。因此，已发布候选的完整评价、比较与复现满足该任务的计算实验范围。GRU-VAE 训练是扩展，原 checkpoint 缺失限制原生成过程重放，但不阻断 GV02-03 的冻结候选评价。

原文没有要求通过湿实验确认功能。当前名单只能解释为计算筛选结果，参考仍为暂定版。家族建模和 Pfam/TM 的发布证据已核验，完整外部建模与扫描尚未重新执行。

真实成员信息、会议、工时及实际答辩记录应来自真实活动；仓库整理不生成这些记录。原过程管理要求保存在 [process-requirements.docx](course/process-requirements.docx)。
