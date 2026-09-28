# 方法与运行入口

## 评价设计

项目先执行 QC、竞争性家族判断及 PF00741 域门槛，再在相同合格池上比较四个目标。Pfam 域同时覆盖部分 GvpJ，因此单独域命中不能代替家族判断。统计保守位点由固定参考 MSA 定义，不由待筛选候选反向决定。

全局比对采用 BLOSUM62、gap open=-10、gap extend=-0.5，并要求两条序列覆盖均 ≥0.80、一致率 ≥0.20、比对分数 ≥0、对齐残基数 ≥20。有效一致率是相同残基数除以较长序列长度，距离为其补数。未达证据条件的距离保留 NaN。

新颖性是到可靠天然参考的最近距离。候选独特性是到 136 条联合合格候选中其他条目的已解析距离均值。排名使用至少 50% 可靠配对比例，25%/50%/75% 三档用于敏感性分析。集合多样性独立计算，并同时报告配对覆盖、条件均值和未知距离上下界。

三种策略均使用同一 114 条排名候选。等权加权、Pareto 非支配层加拥挤距离、分维度轮转均输出固定预算；并列最终由稳定 ID 确定。消融保持硬门槛，只移除一个排序目标。

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
| `vae.py`、`generation.py`、`generation_data.py` | VAE 与数据划分 |

历史模块名和冻结目录名保留，以便与源提交和清单逐项核对。
