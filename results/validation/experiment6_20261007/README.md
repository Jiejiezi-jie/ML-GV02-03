# 实验 6 验证记录

日期：2026-10-07。新增 Pareto 前沿序列特征分析，并更新正式项目报告。原始评价批次与四维评分保持不变。

## 已执行检查

1. 全套测试：333 项通过，零失败、零错误。命令为 `python -m pytest tests -q -W error::pytest.PytestUnhandledThreadExceptionWarning --junitxml=results/validation/experiment6_20261007/pytest_results.xml`；Windows 环境设置 `PYTHONUTF8=1`、`PYTHONIOENCODING=utf-8`。
2. 原评价目录 44 项输出全部通过其原始 SHA-256 清单核验。
3. 新实验在两个独立输出目录运行，24 个产物及清单逐字节一致；见 `reproducibility.json`。
4. 新实验记录的输入及三份实现代码哈希与当前文件一致。114 条候选的 Pareto 层级与主报告一致。
5. 主报告八项产物及生成代码通过清单核验。在临时目录连续生成及替换报告，结果和清单保持一致，预置的无关文件保持不变。
6. 三份新图通过面板对齐、最终 PDF 字体、遮挡检查及 PDF 渲染目检。静态预检的五项非阻断提示见[图表检查说明](../../gv02_03_v2/pareto_analysis_member_a_v1/qa/README.md)。

`summary.json` 记录核验范围，`manifest.json` 保存本目录记录和图表检查记录的哈希。图表工具依赖不属于最小实验环境。

## 范围

本轮重跑新增实验与测试、重新生成报告，并核验原发布评价。2026-09-28 的 QC、全局比对及完整评价复现仍属于当时的验证记录；本轮没有重新训练原 VAE、重放其原始采样，或执行完整家族/Pfam/TM 扫描。没有新增功能实验。已有演示稿保留原内容，实验 6 的新增结果以当前报告为准。
