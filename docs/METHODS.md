# 方法与运行入口

## 评价设计

项目先执行 QC、竞争性家族判断及 PF00741 域门槛，再在相同合格池上比较四个目标。Pfam 域同时覆盖部分 GvpJ，因此单独域命中不能代替家族判断。统计保守位点由固定参考 MSA 定义，不由待筛选候选反向决定。

全局比对采用 BLOSUM62、gap open=-10、gap extend=-0.5，并要求两条序列覆盖均 ≥0.80、一致率 ≥0.20、比对分数 >0、对齐残基数 ≥20。有效一致率是相同残基数除以较长序列长度，距离为其补数。未达证据条件的距离保留 NaN。

新颖性是到可靠天然参考的最近距离。候选独特性是到 136 条联合合格候选中其他条目的已解析距离均值。排名使用至少 50% 可靠配对比例，25%/50%/75% 三档用于敏感性分析。集合多样性独立计算，并同时报告配对覆盖、条件均值和未知距离上下界。

三种策略均使用同一 114 条排名候选。等权加权、Pareto 非支配层加拥挤距离、分维度轮转均输出固定预算；并列最终由稳定 ID 确定。消融保持硬门槛，只移除一个排序目标。

## 选做实验 6：Pareto 前沿与序列特征

该实验已补齐前沿成员的目标取舍解释和序列特征比较，与实验 7 随机筛选对照共同构成两项选做。完整分析方案见 [EXPERIMENT_6_PROTOCOL.md](EXPERIMENT_6_PROTOCOL.md)，结果见 [实验 6 报告](../results/gv02_03_v2/pareto_analysis_member_a_v1/EXPERIMENT_6_REPORT.md)。

### 前沿与目标优势

固定使用原发布的 114 条排名候选及其四维分数，四个目标均最大化。若一个候选在全部目标上不低于另一个，且至少一项更高，则前者支配后者；逐层移除非支配成员，得到 13 层，其中第一前沿有 10 条。前沿按四维共同定义，二维投影不重新定义其成员。

逐条保留原始分数，并以 `100 × 平均并列升序排名 / 114` 计算每一维的同池百分位，最高和最低百分位的所有并列维度均列出。百分位用于解释不同量纲分数的相对位置，不表示功能概率；即使保守性满分，也可能因为满分候选较多而得到较低的平均并列百分位。前沿在目标上的突出表现来自其筛选定义，不算独立有效性验证。

### 预定义特征与对照

主比较为第一前沿 10 条与同池其余 104 条。分析前固定六项序列描述，并完整输出 20 种标准氨基酸的组成比例：

| 特征 | 定义 |
| --- | --- |
| 长度 | 残基总数 |
| 组成 Shannon 熵 | `−Σ p(a) log2 p(a)`，仅对非零比例求和，单位为 bits |
| 最高单一残基比例 | 20 种残基比例的最大值 |
| DEKR 比例 | D、E、K、R 的计数之和除以长度 |
| G/P 比例 | G、P 的计数之和除以长度 |
| 最长同聚物段 | 连续相同残基的最长段长 |

两组分别报告样本数、均值、中位数、四分位数和完整范围；组间报告前沿减对照的均值差及 Cliff's delta。后者为全部跨组配对中 `P(前沿值 > 对照值) − P(前沿值 < 对照值)`，相等贡献为零；不是一对一配对差的统计量。这些特征没有统一的“越高越好”方向，DEKR 比例不表示特定 pH 下的净电荷，G/P 比例也不是结构或功能预测。20 种比例之和为 1，因此各残基差异不是独立效应。

### 长度匹配与等价对照敏感性

敏感性对照按 ID 排序后，用 SciPy `linear_sum_assignment` 为每条前沿候选选取一条不同的非前沿候选，使 10 对的绝对长度差总和最小；不重复使用对照。当前最优总差为 2 个残基。输出每对 ID、长度差、各特征差，以及配对差的汇总。该对照仅减少长度不平衡，不代表随机分组或已消除其他偏差。

长度最优解可能不唯一，因此在首次比较后追加 100 次等价最优匹配检查，固定 `seed=42`。每条边的代价为 `绝对长度差 × (10+1) + Uniform[0,1)`，并在每次重新匹配后断言原始总长度差仍为 2。全部次级随机代价之和小于 10，不能抵消主目标一个残基所对应的 11 个代价单位。所有解和全部特征的均值差均保留，汇总差值范围及负、零、正的次数；零的判断容差为 `1e-12`。

这项事后检查描述对照选择是否影响结论，不是 Bootstrap、独立重复或对全部最优解的均匀抽样，也未穷举全部最优解。差值范围仅适用于这 100 次检查，不作为置信区间；不计算 p 值，不挑选有利的匹配结果。

### 距离证据与复现

逐条保留人工复核标记和 `uniqueness_resolved_fraction`。该覆盖率的分母是原 136 条联合合格池中除自身外的 135 条候选，分子是其中可靠距离的数量；组间报告的是这些逐候选覆盖率的分布。它不同于第一前沿内部的 45 对距离覆盖率，也不同于 Top-20 内部的 190 对覆盖率。未解析距离仍保留缺失状态。

```bash
python experiments/run_pareto_analysis.py --output results/reproductions/pareto_analysis
```

输出目录必须为新的或空目录。入口核验原发布结果清单、分数与 FASTA 的一致性，并检查重算的 Pareto Top-20 与已发布名单相同；新表、图、报告及实现版本由独立哈希清单记录。当前分析只描述单一生成批次中的固定筛选池，候选不构成独立生物学重复；共有组成特征和匹配后的组间差异均不能证明真实功能或因果机制。

## 标准复现

```bash
python experiments/run_full_experiment.py --output results/reproductions/evaluation
python experiments/reproduce_project.py --output results/reproductions/project
```

完整环境与严格比较的使用范围见 [REPRODUCIBILITY.md](REPRODUCIBILITY.md)。

## 单项 QC 和比对

安装扩展依赖后可运行独立入口。默认配置指向当前候选批次；单项检查不自动改变冻结发布结果。

```bash
python experiments/run_quality_audit.py --output-dir results/reproductions/qc
python experiments/run_nearest_reference.py --output-dir results/reproductions/nearest
python experiments/run_qc_report.py --output-dir results/reproductions/qc_report
python experiments/run_similarity.py --config configs/similarity.yaml
```

单项 QC 只根据序列检查模式与长度上限，原始 EOS 结束状态需使用生成元数据核验。完整入口会执行元数据核对。人工预警是复核提示，不自动改写 QC 或家族标签。

自定义候选时复制配置到新文件，修改输入 FASTA、标签与输出目录；保留同批次 ID、原始输入、比对规则和哈希。不同家族应分池比较，新生成池需要重新获得家族、域和保守位点证据，不能直接套用旧的冻结源提交。

## 家族及域证据重建

`configs/family_classification.yaml` 定义参考家族鉴定，`configs/family_classification_vae_member_a_v1.yaml` 定义已发布生成候选的家族鉴定。实现位于 `family_pipeline.py`、`classification.py`、`adjudication.py`。

这些步骤需要 Linux/WSL 中的 HMMER、MAFFT、CD-HIT 和 BLAST+。`environment.yml` 是外部序列工具环境配置；最小冻结复现应使用锁定的 `requirements-*.txt`。完整 Pfam 扫描及 PureseqTM 还需相应数据库和程序，既有结果清单记录原版本。

重建实验前复制配置，并为 processed/model/result 输出使用新的目录，避免覆盖发布证据：

```bash
python experiments/run_family_classification.py --config path/to/new_family_config.yaml
python experiments/run_full_experiment.py --workflow legacy --config path/to/new_scoring_config.yaml
```

第二条命令中的 `legacy` 是保留的命令行兼容名称，表示从序列调用外部工具重新计算原始四维分数；可参考 `configs/gv02_03_vae_member_a_v1.yaml`。发布主流程还会用全局距离替代旧局部比对距离，并严格处理未解析配对。两种结果不能混为同一个研究批次。

## GRU-VAE 生成扩展

固定参考按同源簇隔离，训练/验证/测试为 277/33/36 条。模型为单层 GRU 编码器/解码器，embedding=32、hidden=128、latent=32，使用独立 PAD/BOS/EOS、KL warmup、验证早停及断点恢复。

当前参考暂定，训练入口要求显式声明允许暂定数据。以下命令训练新模型并采样新批次：

```bash
python experiments/train_generator.py --allow-provisional-data --output-dir results/reproductions/vae
python experiments/generate_candidates.py --checkpoint results/reproductions/vae/best.pt --allow-provisional-data --output-dir results/reproductions/generated
```

参数 `--device cpu` 或 `--device cuda:0` 指定设备。代码默认使用当前配置 `generation_member_a_v1.yaml`；采样必须明确传入真实 checkpoint。没有原权重时重新训练不保证得到原来的 1,000 条候选。

## 实现位置

| 模块 | 作用 |
| --- | --- |
| `quality.py`、`sequence_patterns.py` | QC 与序列模式 |
| `classification.py`、`family_pipeline.py` | 家族鉴定 |
| `similarity.py`、`nearest_reference.py` | 全局距离与最近参考 |
| `c_handoff.py`、`c_domains.py`、`c_tm.py`、`c_final_handoff.py` | 批次审计、域/跨膜证据与统一复核 |
| `metrics.py`、`analysis.py`、`selection.py` | 原始评分、统计实验与策略 |
| `d_evaluation.py`、`release.py`、`reproducibility.py` | 冻结评价、验收、数值与哈希比较 |
| `pareto_features.py`、`experiments/run_pareto_analysis.py` | 四维前沿解释、序列特征比较与长度匹配敏感性 |
| `vae.py`、`generation.py`、`generation_data.py` | VAE 与数据划分 |

历史模块名和冻结目录名保留，以便与源提交和清单逐项核对。
