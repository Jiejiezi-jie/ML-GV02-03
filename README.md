# GV02-03：GvpA 候选序列的多目标评价与筛选

本仓库实现课程项目 Track 1 / GV02-03。系统接收已有蛋白质生成候选，使用约束满足度、统计保守性、新颖性和多样性四个目标进行可解释评价，对比等权加权、Pareto 和各维度轮转三种筛选策略，并完成消融、权重/阈值鲁棒性与随机基线实验。下列结果是首次汇报时的旧版基线；根据汇报后的教师意见，项目 V2 将增加序列 VAE 生成、质量控制和竞争性家族鉴定，实施计划见 [改进计划](docs/GV02-03_待做内容与改进计划.md)。

## 教师反馈后的 V2 进度

- 已冻结首选模型为序列 VAE、随机种子 42、生成数 1000 以及 Top-10/20/50 预算；
- 成员 A 交接的 346 条保守暂定 GvpA 参考已按 31 个 70% identity 同源簇划分为
  训练 277 条/16 簇、验证 33 条/9 簇、测试 36 条/6 簇，同源簇无跨集合泄漏；
- 旧版 1224 条/418 簇划分仅作历史开发基线，本轮 VAE 不再用它训练；
- 已实现与冻结配置一致的 GRU 序列 VAE，包括 PAD-aware 编码、KL warm-up、验证/早停、
  checkpoint 严格恢复、逐样本 EOS 生成及完整生成元数据；训练记录包含困惑度、token/EOS
  准确率，生成Manifest包含唯一率、重复、长度和训练/参考完全匹配诊断；
- 已在保守暂定数据上训练至早停并生成 1000 条候选：1000 条全部唯一，972 条自然 EOS；
- QC 硬规则 1000/1000 通过；竞争性家族鉴定得到 147 条主池、47 条 A/J 歧义和 806 条排除；
- 主池中 136 条通过严格域门槛，已完成等权、Pareto、维度轮转的 Top-10/20/50、消融、鲁棒性和随机基线。

阶段说明见 [M0–M1 生成数据准备](reports/M0_M1_生成数据准备.md)和
[成员 B 的 VAE 开发版说明](reports/V2_B_VAE开发版.md)和
[完整 VAE 实验报告](reports/V2_B_VAE完整实验.md)。成员 C 接续 QC 与相似度工作时，
请使用 [成员 B → C 候选序列交接说明](reports/成员B_to_C_候选序列交接说明.md)中冻结的
输入、字段约定、复现命令和验收清单。

## 成员 B 的 VAE 参考代码

组员提供的原始最小参考包保留在 `gv/`，仅作为来源记录。项目接入后的实现位于
[`src/gv_eval/vae.py`](src/gv_eval/vae.py)和
[`src/gv_eval/generation.py`](src/gv_eval/generation.py)，训练与生成入口分别为
[`experiments/train_generator.py`](experiments/train_generator.py)和
[`experiments/generate_candidates.py`](experiments/generate_candidates.py)。

最新交接请先阅读 [成员 A 完整工作记录](reports/成员A_完整工作过程与交接记录.md)
和 [成员 B 实施指南](reports/成员B_VAE实施与数据交接指南.md)。
独立保守开发包位于 `data/processed/gv02_03_v2/member_b_handoff/`，346 条参考已按簇重新划分；
它明确属于暂定开发版本，不能代替仍需核查的最终科学参考。

成员 A 的参考核验和竞争分类入口为 `python -B experiments/run_family_classification.py`，
验收入口为 `python -B experiments/verify_family_deliverables.py`。
完整复现还须依次运行裁决、保守发布包和开发划分脚本，且**顺序不可调换**（验收表会被后续步骤按哈希钉住），
六步顺序见下文方法文档。
交付路径、PF00741 官方边界、规则、结果和未解决的负对照冲突详见
[成员 A 方法与验收](reports/成员A_家族鉴定方法与验收.md)。
当前为计算支持参考，科学验收仍有 14 条名义 GvpJ 对照冲突待核查；
不能把现有生成训练划分视为已经完成高可信参考替换。

训练脚本默认拒绝暂定数据；仅开发验证可显式加入 `--allow-provisional-data`。
模型权重和临时候选不提交 Git，正式 checkpoint 必须等保守参考发布状态明确后重新训练。

## 关键结果

- 候选：T05 生成的 200 条序列；天然参考清洗后 1224 条，CD-HIT 90% 聚类为 418 条代表。
- 资源：官方 PF00741.24 / Gas_vesicle profile；23 个仅由天然参考 MSA 定义的统计保守位点。
- 域证据：37/200 条有可报告 PF00741 命中，27/200 条通过 GA 25 bits + 模型覆盖 0.95 门槛。
- 目标冲突：约束–新颖性 Pearson -0.7961，保守性–新颖性 -0.9127；约束–保守性 0.7498。
- Top-20 等权加权：约束 0.6626、保守性 0.8739、域通过率 100%，但集合多样性 0.6121。
- Top-20 Pareto/轮转：集合多样性 0.8753/0.8857，新颖性 0.5971/0.5893，但域通过率降到 60%/50%。
- 200 次±20%权重扰动：Top-20 Jaccard 均值 0.9064、最小 0.6667；完整排名 Spearman 均值 0.9959。
- 自动化测试：当前全项目 95 项全部通过。

这些都是计算代理结果，不是候选功能准确率或湿实验成功率。PF00741 也并非 GvpA 专一；候选进入下游前仍需结构与实验验证。

## 数据

唯一待筛选输入为 [generated_gvp.fasta](data/candidates/t05/generated_gvp.fasta)，共 200 条。天然参考来自：

- `data/raw/t05/gvpa/GvpA_RefSeq.fasta`
- `data/raw/t05/gvpa/GvpA_NotPartial.fasta`
- `data/raw/t05/gvpa/rescued_GvpA_candidates.fasta`
- `data/raw/design/GvpA.fasta`

参考清洗要求显式 GvpA 标签，排除 GvpJ、partial/fragment/predicted、非标准氨基酸和 50–180 aa 以外记录。混合多家族的 `data/raw/t05/real_gvp.fasta` 不作为纯 GvpA 参考。详情见 [数据说明](data/README.md)和[参考映射](data/processed/gv02_03/reference_mapping.tsv)。

## 环境与完整复现

推荐使用仓库提供的 Conda 配置：

```bash
cd /home/user/wangyuhan/ML-GV02-03
conda env create -f environment.yml
conda activate ml-gv02-03
python -B experiments/run_full_experiment.py --config configs/gv02_03.yaml
python -m pytest -q
```

本次实际环境位于 `/home/user/wangyuhan/envs/ml-gv02-03`；未激活环境时可运行：

```bash
PATH=/home/user/wangyuhan/envs/ml-gv02-03/bin:$PATH \
  /home/user/wangyuhan/envs/ml-gv02-03/bin/python -B \
  experiments/run_full_experiment.py --config configs/gv02_03.yaml
PYTHONPATH=src /home/user/wangyuhan/envs/ml-gv02-03/bin/python -m pytest -q
```

首次运行从 EMBL-EBI InterPro 官方接口下载 PF00741 HMM，下载失败或内容校验失败会停止，不会换成伪造分数。参数和随机种子 42 固定在 [配置文件](configs/gv02_03.yaml)，输入/输出哈希、工具与 profile 版本见 [实验清单](results/gv02_03/manifest.json)。流水线不会修改原始 FASTA。

## 可选候选 QC（新增）

运行 `python -B experiments/run_full_experiment.py --config configs/gv02_03_qc.yaml`
可在评分前启用基础候选质控，输出到 `data/processed/gv02_03_qc/` 和
`results/gv02_03_qc/`。原配置及上面的历史实验结果保持为旧版基线。

QC 硬失败包括空序列、非法氨基酸、长度不在 50–180 aa、完全重复候选
（重复组全部排除），以及低复杂度代理规则：不同残基少于 8 种或某一残基
比例超过 0.35。这些是可配置的初始筛查阈值，尚未经过生物学校准。
`generation_length_cap` 若设置，达到或超过该长度会产生警告而不单独排除；
旧候选缺少生成上限记录，默认 `null` 表示未检查，不能据此认定没有截断。

`candidate_qc.tsv` 保留所有候选的状态和原因，`eligible_candidates.fasta`
保存通过者，评分和 Top-K 表携带 QC 字段。通过数量不足现有分析预算时
明确报错并保留 QC 表，需调整预算后重新运行。原始 FASTA 不被修改。
本次完成改进计划“工作二”的基础序列 QC、训练/参考完全匹配检查及组成偏离
告警；疏水片段、额外结构域和生成 EOS 元数据检查仍待实现。QC 通过不代表 GvpA
家族鉴定通过。完整 QC 实验仍需 HMMER、MAFFT、BLAST 和 CD-HIT 环境运行。

### 训练集与天然参考完全匹配检查

`inputs.qc_training_fasta` 和 `inputs.qc_reference_fasta` 分别指定比对 FASTA。
QC 配置默认使用 V2 暂定训练集 `data/processed/gv02_03_v2/generator_data/train.fasta`
以及天然参考 `data/processed/gv02_03/reference_clean.fasta`。
旧 T05 生成模型的真实训练集未知，所以默认的训练匹配只表示“与指定 V2 训练集重复”，
不能解释为旧模型记忆了训练数据。为新生成批次运行时，应指定实际使用的训练 FASTA。

比较完整氨基酸序列（忽略大小写，不做近似匹配）。输出字段：

- `exact_training_match` / `exact_reference_match`：是否找到完全相同的序列；
- `training_match_ids` / `reference_match_ids`：所有匹配 ID 的有序 JSON 数组；
- `training_match_checked` / `reference_match_checked`：是否实际执行了该项检查。

在 `quality` 下独立设置 `exact_training_match_action` 和
`exact_reference_match_action`，取值为 `warn`（默认，仅告警）或 `exclude`
（将匹配者设为 `qc_pass=false`，不进入评分）。警告同样写入 `qc_reasons`。
比对路径缺省或为 `null` 时标记未检查；此时不能启用 `exclude`。
指定文件不存在、为空、含非法序列或重复 ID 时明确报错，不静默跳过。
同一合法序列对应多个不同 ID 时全部保留。完整运行的摘要包含两类匹配数量，
manifest 记录所用比对 FASTA 的 SHA-256；QC 字段随评分与筛选名单传递。

### 氨基酸组成异常检查

组成检查使用天然参考中每条序列的20种氨基酸频率，计算平均组成，并以
Jensen–Shannon divergence 衡量每条序列相对该平均组成的偏离程度。阈值不是
人工拍脑袋指定的固定距离，而是参考序列自身距离分布的可配置分位数；当前配置为
99%分位数。超过阈值的候选写入 `composition_outlier_warning`，默认只警告，
也可通过 `composition_outlier_action: exclude` 设置为硬失败。

当前基线来自尚待 GvpA/GvpJ 竞争性核验的名义 GvpA 参考，因此这项结果只能作为
异常组成筛查，不是家族或功能判断。完成最终高可信参考集后必须重新校准阈值。

本配置已在旧版200条候选上完整运行：170条通过基础硬QC，30条因长度超过
180 aa被排除；训练集和天然参考完全匹配均为0条。145条触发组成异常警告，
但不会因此被排除。过滤后的170条中有37条PF00741命中、27条通过当前域门槛。
正式产物位于 `data/processed/gv02_03_qc/` 和 `results/gv02_03_qc/`，具体解释和
复现边界见 [候选质量控制增量报告](reports/M3_候选质量控制增量.md)。

## C 模块：独立全局相似度与距离

运行 `python -B experiments/run_similarity.py --config configs/similarity.yaml`，
输出成对比对证据、双向覆盖率、距离矩阵和对应 ID 顺序。
未通过覆盖率、identity、得分或最小比对长度门槛的序列对标记 `unresolved`，
矩阵中保留 NaN，不赋最大距离。原版 BLAST 评分与结果保持不变。
默认输入为旧 T05 混合候选，仅作开发验证；正式主分析池与歧义池需分别输入。
方法、公式、可靠性门槛、缺失值处理及 D 的接入说明见
[全局相似度模块说明](docs/similarity.md)。本阶段尚未实现快速检索和自动家族分池。

## C：序列模式告警（QC-only）

新增可配置的连续疏水片段、单残基重复、短 motif 串联重复告警，记录片段位置；
`qc_warnings` 独立列出警告，旧 `qc_reasons` 保持兼容。新检测不新增淘汰条件。

```bash
python experiments/run_quality_audit.py --config configs/gv02_03_pattern_qc.yaml --output-dir results/my_qc_audit
```

输出目录须为空，运行不依赖外部生信工具。旧实验配置与结果不覆盖。
阈值尚未正式校准，疏水告警不是跨膜预测；范围与字段见 [QC 模式告警说明](docs/qc-patterns.md)。

## C：最近可靠参考匹配

完整比较候选与指定参考库，在通过全局比对可靠性门槛的参考中选择最小距离，
输出最近参考 ID、全部并列 ID、identity、双向覆盖率及距离；无可靠匹配保留 `unresolved`。

```bash
python experiments/run_nearest_reference.py --config configs/nearest_reference.yaml --output-dir results/my_nearest_reference
```

不改变旧 QC 或评分；默认输入仍是未正式确认家族的开发数据。
字段、并列规则及边界见 [最近参考匹配说明](docs/nearest-reference.md)。

## C：统一 QC 汇总与人工复核清单

合并已生成的 QC 和最近参考证据，生成逐条总表、复核清单、中文 Markdown 报告及审计记录。
输入先核验哈希、ID 和长度；保留原有 `qc_pass`，不新增淘汰条件。

```bash
python experiments/run_qc_report.py --config configs/qc_report.yaml --output-dir results/my_qc_report
```

需要先完成上面的 QC-only 与最近参考运行，详情见 [统一 QC 汇总说明](docs/qc-report.md)。

## C：B 真实 VAE 批次交接

对 B 的冻结批次（1,000 条候选）进行导出元数据/EOS 审计、基础 QC 复算与扩展警告、
346 条参考的全局最近匹配，以及主池 147 / 模糊池 47 的独立距离矩阵。
读取外部 B 快照，不改变原家族分池。原交接运行未纳入额外结构域和跨膜拓扑；
后续已分别补充 [本地 Pfam 域架构复核](results/gv02_03_c_domains_member_a_v1/README.md)
和 [本地跨膜预警检查](results/gv02_03_c_tm_member_a_v1/README.md)，均不改变 B 分池。
运行方式与文件含义见 [真实批次 C 交接说明](docs/c-handoff.md)。
三批证据的逐条汇总及下游边界见 [C→下游统一复核交接](docs/c-final-handoff.md)。

## 四个指标（旧版基线）

| 指标 | 实现 |
| --- | --- |
| 约束满足度 | PF00741 domain bit score 相对参考命中的经验百分位 × 模型覆盖因子；另报 GA+覆盖门槛通过状态 |
| 保守性 | 参考 MSA 中 gap≤0.10、共识≥0.90 的 23 个列上，候选与共识一致的比例；缺口计 0 |
| 新颖性 | `1 - max_reference(pident × alignment_length / max(lengths))` |
| 独特性/集合多样性 | 候选到其他候选的平均距离 / 筛选子集平均两两距离 |

完整逐条证据在 [candidate_scores.tsv](results/gv02_03/tables/candidate_scores.tsv)。

## 目录

```text
ML-GV02-03/
├── configs/gv02_03.yaml             # 冻结参数
├── data/
│   ├── candidates/                   # 原始候选
│   ├── raw/                          # 原始天然来源（不覆盖）
│   └── processed/gv02_03/            # 清洗参考、MSA、位点和距离矩阵
├── models/pfam/                      # PF00741.24 profile 与来源元数据
├── src/gv_eval/                      # I/O、工具封装、指标、策略、分析、流水线
├── experiments/run_full_experiment.py
├── results/gv02_03/
│   ├── raw/                          # HMMER/BLAST 原始结果
│   ├── tables/                       # 评分、相关、策略、消融、鲁棒性、随机基线
│   ├── selections/                   # 三策略 K=10/20/50 的 TSV 与 FASTA
│   ├── figures/                      # 论文图件
│   ├── summary.json
│   └── manifest.json
├── reports/                          # M1–M5、最终 PDF 与答辩 PPT
├── tests/                            # 自动化测试
└── docs/project_management/          # 真实过程材料填写说明与模板
```

`results/gv02_03/work/blastdb/` 和工具日志属于可重建中间文件，不纳入版本控制；正式评分、原始比对表、图件和名单均保留。

## 报告与答辩材料

- [M1 功能需求分析与任务建模](reports/M1_功能需求分析与任务建模报告.md)
- [M2 探索性数据分析](reports/M2_探索性数据分析.md)
- [M3 评价体系与初步筛选](reports/M3_评价体系与初步筛选.md)
- [M4 完整实验与鲁棒性分析](reports/M4_完整实验与鲁棒性分析.md)
- [M5 最终报告](reports/M5_最终报告.md)及 PDF
- `reports/GV02-03_答辩PPT.pptx`

过程评分还要求真实会议、分工、Git 贡献和照片。仓库只提供模板，不会代替成员编造这些证据；项目组应从现在起按 [过程材料说明](docs/project_management/README.md)持续补充。

## 主要限制

- 大多数 T05 候选没有 PF00741 命中，说明输入是混合 Gvp 生成池，而非已确认的纯 GvpA 池。
- 统计保守位点不是实验验证关键位点；PF00741 命中也不能区分全部 GvpA/GvpJ 情形。
- BLAST 局部比对经过覆盖校正，但不是严格全局结构相似性。
- 无功能真值，不能通过调参声称找到“最优蛋白”；推荐将等权质量组与 Pareto/轮转探索组交给独立结构和实验验证。
