# 教师参考材料

本目录只保留 T05 的输入来源和生成方法证据，服务于 **GV02-03：GvpA 候选的四目标评价与筛选**。本项目评价域约束满足度、关键位点保守性、新颖性和多样性，比较筛选策略并检验稳定性，不要求训练新模型。教师指南位于 [project-pool(1).pdf](../project-pool(1).pdf)。

## 保留内容

| 文件 | 用途 |
| --- | --- |
| [t05/T05-report.docx](t05/T05-report.docx) | 原始完整报告，用于追溯数据背景、生成和原评价流程 |
| [t05/doc/dataset.md](t05/doc/dataset.md) | 数据收集与处理说明 |
| [t05/doc/model_training.md](t05/doc/model_training.md) | 原模型训练、生成设置及交付说明 |
| [t05/src/model_training/model.py](t05/src/model_training/model.py) | 已有 Transformer 架构 |
| [t05/src/model_training/generate.py](t05/src/model_training/generate.py) | 原采样入口与参数 |
| [t05/src/model_training/train.py](t05/src/model_training/train.py) | 训练来源证据，用于理解输入和权重产生方式 |
| [t05/src/model_training/gvp_ckpt/model.pt](t05/src/model_training/gvp_ckpt/model.pt) | 原有权重，供后续候选不足时评估能否补充采样 |

原材料内容未修改，目录示例和命令可能引用已删除工程或旧绝对路径。实际入口以 [数据说明](../data/README.md) 为准：候选是 `data/candidates/t05/generated_gvp.fasta`，天然来源及混合家族背景位于 `data/raw/`。

## 生成资源的限制

训练与生成代码缺少 `utils.py`，训练还缺 `dataset.py`，并包含 AutoDL 绝对路径。缺失内容包括精确词表顺序、BOS/EOS/PAD 索引、解码和采样函数。仅知道词表含 23 个 token 不能恢复权重对应的字母映射，不能任意补词表就声称复现。保留 checkpoint 不表示当前可直接运行；本轮没有加载权重、重新训练或补充生成。

源码 `model.py` 默认是 **4 层、4 个注意力头、256 维表示**；生成设置为 temperature=0.8、top-k=10、最大长度 512、批量 200 条。原报告中的“4–6 层、8 个头”等描述与源码不一致，应区分报告描述和实际保留实现。

T05 使用混合 Gvp 数据，200 条生成序列没有亚家族标签。旧报告的统计一致性和结构合理性结论不能直接作为本项目的 GvpA 身份结论或功能标签；“90% 去重”说明也不能替代对数据和相似度算法的实际核验。后续须建立可追溯的参考集、比对、四目标指标与筛选实验。

## 范围修正与已删除内容

原任务记录中的 PF01132 对应 EFP（EF-P 的 OB 域），不能作为 GvpA 域标识。相关的 PF00741／Gas_vesicle 包含 GvpJ，命中不足以单独确认 GvpA 亚家族。原文件保留不改；本项目报告和后续实现采用核验后的标识，并记录数据库版本和亚家族判定规则。[NCBI PF01132](https://www.ncbi.nlm.nih.gov/Structure/cdd/pfam01132)、[NCBI PF00741](https://www.ncbi.nlm.nih.gov/Structure/cdd/pfam00741)

第二轮已删除整个 `references/design/`、design 的 10 条候选和旧评分表、T05 重复 PPT、原 README 及旧预处理/统计分析材料。design 旧生成与五维评分工程含 mock 或启发式实现，不能替代正式评价。具有独有序列的 design 天然 FASTA 继续保留在 `data/raw/design/GvpA.fasta`。

[preprocessing/analyze_inputs.py](../preprocessing/analyze_inputs.py) 负责当前输入审计，结果位于 [results/input_audit/](../results/input_audit/)。`docs/cleanup/` 第一轮记录是历史快照，第二轮清单和 [整理说明](../docs/整理说明.md) 记录本轮精简去向。
