# 结果目录

## 直接使用的最终结果

- [`gv02_03_v2/d_evaluation_member_a_v1/`](gv02_03_v2/d_evaluation_member_a_v1/)：1,000 条审计记录、114 条排名候选、三策略 Top-10/20/50、相关性、消融、鲁棒性及随机对照。
- [`gv02_03_v2/pareto_analysis_member_a_v1/`](gv02_03_v2/pareto_analysis_member_a_v1/)：选做实验 6，四维前沿、10 对 104 的序列特征对照、长度匹配及 100 次等价最优匹配敏感性。入口为 [EXPERIMENT_6_REPORT.md](gv02_03_v2/pareto_analysis_member_a_v1/EXPERIMENT_6_REPORT.md)。
- [`validation/`](validation/)：2026-09-28 交付整理后的测试、完整复现、独立副本验收及 PPT 校验记录，入口为 `summary.json`。
- 正式解读及图表位于 [`../reports/`](../reports/)，阅读入口为 [PROJECT_REPORT.md](../reports/PROJECT_REPORT.md)。

四项必做及选做实验 6、7 均已完成；实验 5 下游性质预测未做。两个正式实验目录使用独立清单，新增特征分析不改变原来的 44 个评价输出。

### 实验 6 文件索引

| 文件 | 内容 |
| --- | --- |
| `candidate_features.tsv`、`front_details.tsv` | 完整 114 条候选特征与 10 条前沿的四维优势、哈希、覆盖及复核标记 |
| `group_comparison.tsv`、`layer_summary.tsv` | 两组比较的全部特征及 13 个 Pareto 层的描述统计 |
| `matched_pairs.tsv`、`matched_feature_differences.tsv` | 一对一长度匹配及逐对特征差 |
| `tie_match_pairs.tsv`、`tie_match_differences.tsv`、`tie_match_summary.tsv` | 100 次等价最优匹配、完整差值及方向稳定性 |
| `summary.json`、`manifest.json`、`figures/` | 实际分析参数、证据边界、文件哈希与可视化 |

主比较发现前沿平均较长；长度匹配后的组成熵差在所检查的 100 个等价解中均为正，最高单一残基比例差均为负，G/P 比例差则方向不一致。这些为当前批次的描述性结果，不作功能或因果推断。

## 当前批次依赖的证据

| 目录 | 内容 |
| --- | --- |
| `gv02_03_v2/family/` | 天然参考家族校准、控制序列及验证摘要 |
| `gv02_03_v2/vae_member_a_v1/` | 当前生成候选的家族判定与初始评分证据 |
| `gv02_03_c_handoff_member_a_v1/` | QC、到天然参考的全局距离、同池距离矩阵及质量预警 |
| `gv02_03_c_domains_member_a_v1/` | Pfam 结构域复核证据 |
| `gv02_03_c_tm_member_a_v1/` | 跨膜扫描复核证据 |
| `gv02_03_c_final_handoff_member_a_v1/` | 汇总域、跨膜与序列质量的统一审计 |
| `gv02_03_v2/integration_20260928/` | 清理前集成快照的复现和 304 项测试记录 |

证据目录的原名称用于关联冻结源提交与哈希清单，因此保留不变。最终筛选分数以 `d_evaluation_member_a_v1/ranking_pool_scores.tsv` 为准。2026-09-28 整理后为 302 项测试，清理时移除了两个旧批次快照测试；2026-10-07 新增前沿分析后，333 项测试通过。日期较早的复现记录保留原始范围，不代表新增实验或外部扫描在该次验证中执行过。

重新生成前沿分析：

```bash
python experiments/run_pareto_analysis.py --output results/reproductions/pareto_analysis
```

本地重新运行的输出统一放在 `reproductions/`，由 Git 忽略。已发布结果由 `manifest.json` 记录文件哈希，请使用新的输出目录重新运行。
