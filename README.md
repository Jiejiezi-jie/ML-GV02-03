# GvpA 多目标评价与筛选

对气囊蛋白 GvpA 的候选氨基酸序列进行质量审计、四维评价和多目标筛选。项目在相同候选池与筛选预算下比较加权评分、Pareto 非支配排序和分维度轮转，输出候选名单、实验结果及可复现报告。

**阅读入口：** [从零理解项目](docs/PROJECT_GUIDE.md) · [实验报告](reports/PROJECT_REPORT.md) · [中期材料](reports/midterm/README.md) · [方法定义](docs/METHODS.md) · [复现说明](docs/REPRODUCIBILITY.md)

## 项目能力

- 审计候选来源、序列质量、家族归属、结构域支持及人工复核标记。
- 计算约束满足度、保守性、新颖性和候选独特性，并保留未解析的距离。
- 在固定预算下生成三种策略的 Top-10、Top-20、Top-50 名单及 FASTA。
- 分析指标关系、策略差异、逐维消融、权重与门槛敏感性。
- 比较同池随机筛选，解释 Pareto 前沿的序列特征及长度匹配敏感性。
- 校验冻结输入与输出哈希，保留版本、配置和逐候选证据。

```text
参考序列与生成候选
       ↓
质量控制、家族判定、结构域证据
       ↓
全局序列比对与四维评分
       ↓
证据充分性检查 → 三策略筛选
       ↓
相关性、消融、敏感性、随机对照、前沿特征分析
       ↓
候选名单、审计表、图表与报告
```

## 快速运行

使用 **Python 3.12 和 Git**。评价与前沿分析使用 CPU，不需要 GPU 或生成模型权重。请完整克隆仓库，冻结入口需要读取 Git 历史中的原始输入。

```bash
git clone https://github.com/Jiejiezi-jie/ML-GV02-03.git
cd ML-GV02-03
python -m venv .venv
```

激活环境：

```powershell
# Windows PowerShell
.venv\Scripts\Activate.ps1
```

```bash
# Linux / macOS
source .venv/bin/activate
```

安装依赖并运行两个实验入口：

```bash
python -m pip install -r requirements-frozen-v2.txt
python experiments/run_full_experiment.py --output results/reproductions/evaluation
python experiments/run_pareto_analysis.py --output results/reproductions/pareto_analysis
```

- 第一个入口核验冻结输入，重算评价、筛选、消融、鲁棒性和随机对照，并独立重复检查 44 个输出。
- 第二个入口从已发布的同批次评价结果重算前沿、特征及匹配分析，生成 24 个输出。
- 每次运行使用新的或空目录；产物默认不覆盖已发布结果。下载 ZIP 或缺少源提交的浅克隆不能代替完整克隆。

生成报告：

```bash
python experiments/build_project_report.py --output results/reproductions/report
```

报告生成入口读取并核验已发布评价和前沿分析结果。QC、全局比对的扩展复现、测试与自定义输入详见[复现说明](docs/REPRODUCIBILITY.md)。

## 发布结果

发布批次采用固定种子 42，保守天然参考包含 346 条序列。

| 阶段 | 候选数 |
| --- | ---: |
| 生成 / 基础质量检查通过 | 1,000 / 1,000 |
| GvpA 家族支持 | 147 |
| 质量、家族、结构域联合门槛通过 | 136 |
| 四维证据充分，进入排名 | 114 |
| 第一 Pareto 前沿 | 10 |

三种策略在同一 114 条候选池上的 Top-20 结果：

| 策略 | 约束满足度 | 保守性 | 新颖性 | 候选独特性 |
| --- | ---: | ---: | ---: | ---: |
| 等权加权 | 0.8967 | 0.9891 | 0.3073 | 0.4051 |
| Pareto | 0.6070 | 0.9739 | 0.3751 | 0.4601 |
| 分维度轮转 | 0.6254 | 0.9761 | 0.3121 | 0.4179 |

加权筛选更重视约束支持，Pareto 保留更高的新颖性和候选独特性。200 次权重扰动的名单 Jaccard 均值为 **0.9385**；27 组门槛组合的范围为 **0.6000–1.0000**。

前沿特征分析显示：前沿候选平均更长；长度匹配后，组成熵更高、单一残基最高占比更低的方向在检查的 100 组等价最优匹配中一致，G/P 比例差的方向不一致。完整结果及解释见[项目报告](reports/PROJECT_REPORT.md)和[前沿分析报告](results/gv02_03_v2/pareto_analysis_member_a_v1/EXPERIMENT_6_REPORT.md)。

![三种策略的质量与集合多样性](reports/figures/strategy_comparison.png)

## 交付文件

| 文件或目录 | 用途 |
| --- | --- |
| [docs/PROJECT_GUIDE.md](docs/PROJECT_GUIDE.md) | 从基础概念到逐项实验的自学教材，解释方法、原因、结果和推论 |
| [reports/PROJECT_REPORT.md](reports/PROJECT_REPORT.md) | 正式研究问题、结果、适用场景与局限 |
| [reports/midterm/](reports/midterm/README.md) | 中期答辩 PPT、讲稿及已有项目概览 |
| [docs/METHODS.md](docs/METHODS.md) | 指标、比对规则、筛选算法与实现位置 |
| [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) | 环境、运行命令、验证范围与历史版本 |
| [results/README.md](results/README.md) | 候选名单、数据表、图表及验证记录索引 |
| [docs/COURSE_REQUIREMENTS.md](docs/COURSE_REQUIREMENTS.md) | GV02-03 原始要求与交付内容对应关系 |

```text
configs/        实验参数及工具配置
data/           参考序列、候选、划分和输入证据
preprocessing/  输入审计
models/         家族/Pfam 模型与生成记录
src/gv_eval/    评价、比对、筛选及生成扩展实现
experiments/    实验、复现及报告命令行入口
results/        正式实验结果与验证记录
reports/        正式实验报告、图表与中期材料
docs/           自学教材、方法、复现说明及课程原件
references/     数据来源及家族争议复核证据
tests/          算法、数据完整性与端到端测试
```

## 适用范围

本项目完成 GV02-03 的四项必做实验，并选择实验 6（Pareto 前沿分析）与实验 7（随机筛选对照）作为两项选做。实验 5 的下游性质预测未实现。结题答辩 PPT 尚未制作，已有中期材料按原版本保留。

结果用于候选排序和后续验证规划。参考集仍为暂定版本，尚无表达、结构、组装或功能实验确认。未知距离保留为缺失值；集合多样性同时报告配对覆盖率及上下界，不能仅凭已解析对均值判断整体优劣。

GRU-VAE 训练与生成代码作为扩展保留。原始权重未交付，原训练与采样尚不能重放；冻结候选的评价和筛选不依赖该权重。完整家族建模、Pfam 和 PureseqTM 重扫需要额外工具及数据库，现有复现流程会明确区分重新计算与已有证据核验。
