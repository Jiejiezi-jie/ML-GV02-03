# GvpA 多目标评价与筛选

面向气囊蛋白 GvpA 候选序列的计算评价工具。项目结合序列质量控制、竞争性家族鉴定、结构域证据和全局序列比对，在约束满足度、保守性、新颖性与多样性之间比较筛选策略，输出可追溯的候选名单。

[项目讲解](docs/PROJECT_GUIDE.md) · [项目报告](reports/PROJECT_REPORT.md) · [演示文稿](reports/PROJECT_PRESENTATION.pptx) · [复现说明](docs/REPRODUCIBILITY.md) · [方法说明](docs/METHODS.md) · [课程要求对应表](docs/COURSE_REQUIREMENTS.md)

中期汇报材料：[汇报 PPT](reports/midterm/机器学习中期汇报.pptx) · [10 分钟逐页演讲稿](reports/midterm/GV02-03_逐页演讲稿_10分钟.docx)

## 功能

- **候选审计**：保留每条序列的来源、SHA-256、质量检查、家族归属及人工复核标记。
- **四维评价**：计算结构域支持、统计保守位点一致性、天然参考距离与候选独特性。
- **多目标筛选**：比较等权加权、Pareto 非支配排序和分维度轮转，导出 Top-10/20/50 的 TSV 与 FASTA。
- **实验分析**：相关性与散点矩阵、Pareto 前沿、逐维消融、权重及门槛敏感性、同池随机基线。
- **可复现运行**：固定输入版本，校验文件哈希，独立复跑比较输出，并保留软件版本和审计清单。
- **候选生成扩展**：提供 GRU 序列 VAE 的训练、断点恢复、采样与元数据记录代码。

## 流程

```text
候选序列与天然参考
        ↓
质量控制、GvpA/GvpJ 家族判定、Pfam 域门槛
        ↓
全局比对、保守位点与四维评分
        ↓
证据完整性检查 → 三策略筛选 → 鲁棒性与随机对照
        ↓
候选名单、逐条审计表、图表与项目报告
```

未解析的距离保留为缺失值。候选独特性用于逐条排序；集合多样性同时报告已解析对均值、配对覆盖率和未知距离对应的上下界。

## 发布结果

本仓库发布一个固定种子为 42 的序列 VAE 批次。

| 阶段 | 序列数 |
| --- | ---: |
| 保守天然参考 | 346 |
| 生成候选 / 基础 QC 通过 | 1,000 / 1,000 |
| GvpA 家族支持 | 147 |
| QC、家族、结构域联合门槛通过 | 136 |
| 四维证据充分、进入排名 | 114 |
| 第一 Pareto 前沿 | 10 |

三种策略均在同一 114 条候选池比较，每种提供 Top-10、Top-20、Top-50。主预算 Top-20 的结果如下。

| 策略 | 约束满足度 | 保守性 | 新颖性 | 候选独特性 |
| --- | ---: | ---: | ---: | ---: |
| 等权加权 | 0.8967 | 0.9891 | 0.3073 | 0.4051 |
| Pareto | 0.6070 | 0.9739 | 0.3751 | 0.4601 |
| 分维度轮转 | 0.6254 | 0.9761 | 0.3121 | 0.4179 |

等权加权更侧重约束支持，Pareto 保留更高的新颖性和候选独特性。200 次 ±20% 权重扰动的 Top-20 Jaccard 均值为 **0.9385**；27 组门槛组合的 Jaccard 范围为 **0.6000–1.0000**。完整分析与 1,000 次随机对照见[项目报告](reports/PROJECT_REPORT.md)。

![三种策略的质量与集合多样性](reports/figures/strategy_comparison.png)

## 快速开始

需要 **Git 和 Python 3.11/3.12**，推荐 Python 3.12。冻结复现入口从历史提交读取已校验的输入，请完整克隆仓库，保留 Git 历史。

```bash
git clone https://github.com/Jiejiezi-jie/ML-GV02-03.git
cd ML-GV02-03
python -m venv .venv
```

激活虚拟环境：

```powershell
# Windows PowerShell
.venv\Scripts\Activate.ps1
```

```bash
# Linux / macOS
source .venv/bin/activate
```

安装依赖并运行：

```bash
python -m pip install -r requirements-frozen-v2.txt
python experiments/run_full_experiment.py --output results/reproductions/evaluation
```

该入口核验冻结输入，重算评分、筛选与全部实验，再独立重复运行，比较 44 个输出产物。评价复现使用 CPU，不需要生成模型权重。每次运行请使用新的或空的输出目录。

### 扩展复现与测试

```bash
python -m pip install -r requirements-integration.txt
python experiments/reproduce_project.py --output results/reproductions/project
```

扩展入口还会两次重算质量控制和全局比对，核验参考、家族、域与跨膜证据，并与发布结果比较。完整测试包含 VAE，需要额外安装 PyTorch：

```bash
python -m pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cpu
python -m pytest tests -q
python experiments/check_frozen_checkout.py --report results/reproductions/checkout.json
```

重新生成报告和图表：

```bash
python experiments/build_project_report.py --output results/reproductions/report
```

自定义输入、训练采样和外部工具配置见[方法说明](docs/METHODS.md)。

## 目录

```text
configs/        实验、生成与序列检查配置
data/           原始参考、生成候选、数据划分与逐条证据
models/         Pfam/家族模型、VAE 训练记录
preprocessing/  FASTA 输入审计
src/gv_eval/    评价、比对、生成与筛选实现
experiments/    命令行入口、复现与报告生成
results/        冻结实验结果及验证记录
reports/        当前批次的项目报告、图表和演示材料
docs/           方法、复现及原始课程要求
references/     数据来源与家族争议复核证据
tests/          单元测试、数据完整性和端到端测试
```

主要结果位于 [`results/gv02_03_v2/d_evaluation_member_a_v1/`](results/gv02_03_v2/d_evaluation_member_a_v1/)，各证据目录的用途见[结果索引](results/README.md)：

| 文件 | 内容 |
| --- | --- |
| `candidate_audit.tsv` | 全部 1,000 条候选及排除原因 |
| `ranking_pool_scores.tsv` / `ranking_pool.fasta` | 114 条排名候选及四维分数 |
| `weighted_sum_top20.*`、`pareto_top20.*`、`dimension_round_robin_top20.*` | 三策略的 Top-20 表格与序列 |
| `strategy_summary.tsv`、`strategy_overlap.tsv` | 策略质量、覆盖率与名单重叠 |
| `ablation_summary.tsv`、`weight_robustness.tsv`、`threshold_coverage_sensitivity.tsv` | 消融及敏感性实验 |
| `random_comparison.tsv`、`manifest.json` | 随机对照与文件哈希 |

## 结果适用范围

参考集仍为保守筛选后的暂定版本；计算名单用于确定后续验证优先级，尚无表达、结构、组装或功能实验确认。原始 VAE 权重未随批次交付，因此原训练与原采样尚未重放；冻结候选的评价与筛选可完整复现。家族建模和完整 Pfam/PureseqTM 扫描的发布证据已核验，重新执行这些步骤需要额外工具与数据库。

本项目对应《机器学习综合实践》GV02-03 多目标候选评价。原始任务的四项必做实验、研究问题及交付文件对应关系见[要求核对表](docs/COURSE_REQUIREMENTS.md)。
