# 结果目录

## 直接使用的最终结果

- [`gv02_03_v2/d_evaluation_member_a_v1/`](gv02_03_v2/d_evaluation_member_a_v1/)：1,000 条审计记录、114 条排名候选、三策略 Top-10/20/50、相关性、消融、鲁棒性及随机对照。
- [`validation/`](validation/)：本次交付整理后的测试、完整复现、独立副本验收及 PPT 校验记录，入口为 `summary.json`。
- 正式解读及图表位于 [`../reports/`](../reports/)，阅读入口为 [PROJECT_REPORT.md](../reports/PROJECT_REPORT.md)。

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

证据目录的原名称用于关联冻结源提交与哈希清单，因此保留不变。最终筛选分数以 `d_evaluation_member_a_v1/ranking_pool_scores.tsv` 为准。当前测试为 302 项，清理时移除了两个只针对旧批次文件快照的测试，算法测试继续保留。

本地重新运行的输出统一放在 `reproductions/`，由 Git 忽略。已发布结果由 `manifest.json` 记录文件哈希，请使用新的输出目录重新运行。
