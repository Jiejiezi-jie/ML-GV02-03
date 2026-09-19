# GV02-03：GvpA 候选序列的多目标评价与筛选

本仓库实现课程项目 Track 1 / GV02-03。系统接收已有蛋白质生成候选，使用约束满足度、统计保守性、新颖性和多样性四个目标进行可解释评价，对比等权加权、Pareto 和各维度轮转三种筛选策略，并完成消融、权重/阈值鲁棒性与随机基线实验。下列结果是首次汇报时的旧版基线；根据汇报后的教师意见，项目 V2 将增加序列 VAE 生成、质量控制和竞争性家族鉴定，实施计划见 [改进计划](docs/GV02-03_待做内容与改进计划.md)。

## 教师反馈后的 V2 进度

- 已冻结首选模型为序列 VAE、随机种子 42、计划生成 1000 条候选以及 Top-10/20/50 预算；
- 已将 1224 条暂定 GvpA 参考按 418 个 CD-HIT 90% 同源簇划分为训练 978 条、验证 123 条、测试 123 条，同源簇之间无跨集合泄漏；
- 已保存固定词表、逐序列来源与哈希及划分清单，运行方式为 `python -B experiments/prepare_generator_data.py`；
- 尚未完成 GvpA/GvpJ 竞争性核验，因此该划分明确标记为暂定数据，不能作为最终高可信训练集，也尚未开始正式 VAE 训练。

阶段说明见 [M0–M1 生成数据准备](reports/M0_M1_生成数据准备.md)。

## 关键结果

- 候选：T05 生成的 200 条序列；天然参考清洗后 1224 条，CD-HIT 90% 聚类为 418 条代表。
- 资源：官方 PF00741.24 / Gas_vesicle profile；23 个仅由天然参考 MSA 定义的统计保守位点。
- 域证据：37/200 条有可报告 PF00741 命中，27/200 条通过 GA 25 bits + 模型覆盖 0.95 门槛。
- 目标冲突：约束–新颖性 Pearson -0.7961，保守性–新颖性 -0.9127；约束–保守性 0.7498。
- Top-20 等权加权：约束 0.6626、保守性 0.8739、域通过率 100%，但集合多样性 0.6121。
- Top-20 Pareto/轮转：集合多样性 0.8753/0.8857，新颖性 0.5971/0.5893，但域通过率降到 60%/50%。
- 200 次±20%权重扰动：Top-20 Jaccard 均值 0.9064、最小 0.6667；完整排名 Spearman 均值 0.9959。
- 自动化测试：16 项全部通过。

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
本次仅完成改进计划“工作二”的基础序列 QC；训练/参考完全匹配、组成偏离、
疏水片段、额外结构域和生成 EOS 元数据检查仍待实现。QC 通过不代表 GvpA
家族鉴定通过。完整 QC 实验仍需 HMMER、MAFFT、BLAST 和 CD-HIT 环境运行。

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
