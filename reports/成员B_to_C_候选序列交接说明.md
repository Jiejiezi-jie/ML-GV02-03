# 成员 B → 成员 C：VAE 候选序列交接说明

- 日期：2026-09-23
- 交接批次：`member_a_v1_seed42`
- 代码分支：`feat/sequence-vae`
- 成员 B 基线提交：`4e85ed4`

## 1. 交接结论

成员 B 已使用成员 A 提供的暂定保守 GvpA 数据训练 GRU 序列 VAE，并冻结一批
1000 条候选序列及逐条生成元数据。成员 C 可以直接从本文件列出的 FASTA、TSV 和
Manifest 开始质量控制，不需要 checkpoint，也不需要重新训练或重新生成。

本批次状态为 `development_provisional`，上游参考 release 状态为
`provisional_not_scientifically_final`。因此它适合完成 QC、相似度模块开发和可复现实验，
不能表述为已获得生物学功能确认的最终 GvpA 候选。

## 2. 冻结输入

以下三个文件共同构成本次 B → C 交接，不可只复制 FASTA：

| 文件 | 用途 | SHA-256 |
| --- | --- | --- |
| `data/generated/sequence_vae/member_a_v1_seed42/vae_candidates.fasta` | 1000 条候选序列 | `d93b92170977fad475104f58c3a6579d2dd6bee4d99c935d94978db454af7d88` |
| `data/generated/sequence_vae/member_a_v1_seed42/generation_metadata.tsv` | 逐条生成来源、停止原因、长度和序列哈希 | `60e95e6647b6e417bc32215f1e104e8eb37a5eba1c7692c9a8a0d143979e17a2` |
| `data/generated/sequence_vae/member_a_v1_seed42/generation_manifest.json` | 配置、输入输出哈希、环境和汇总统计 | `415fd5a1339440ea1cf94de49faa3ca55680e82fec3d3ea16389e23d04876c72` |

关键来源标识：

- 随机种子：42；
- 参考 release：`member_a_v1_conservative_provisional`；
- split manifest SHA-256：
  `28e2c439d1946a789cb3a3812d4c54f1cce0502dbfa9429d77c238425b501f30`；
- 最佳 checkpoint SHA-256：
  `d2801a2a130928739f44f9726bcb5bbd7be229fadfe2f05fba5d0f6cc87debdf`。

checkpoint 本身按仓库策略不提交 Git，但 C 的工作不依赖 checkpoint。候选 FASTA、元数据和
Manifest 已提交，足以对本批次进行逐条审计。

## 3. 元数据字段

`generation_metadata.tsv` 与 FASTA 使用完全相同的 `sequence_id` 和顺序。C 合并数据时必须以
`sequence_id` 为键，并同时核对 `sequence_sha256`，不得依赖行号或手工复制。

| 字段组 | 字段 | 含义 |
| --- | --- | --- |
| 身份 | `sequence_id`, `sample_index`, `latent_id` | 候选稳定 ID、生成次序和潜变量样本编号 |
| 模型 | `generator_type`, `model_checkpoint`, `checkpoint_sha256` | 生成器类型及 checkpoint 来源；路径仅作记录 |
| 数据来源 | `reference_release_id`, `split_manifest_sha256`, `vocabulary_sha256` | 上游数据和词表版本 |
| 采样 | `generation_seed`, `temperature`, `top_k`, `top_p` | 固定生成参数 |
| 停止状态 | `raw_token_length`, `terminated_by_eos`, `hit_generation_cap`, `stop_reason` | EOS 或长度上限停止证据 |
| 序列审计 | `sequence_length`, `sequence_sha256` | FASTA 长度和内容哈希 |

`terminated_by_eos` 与 `hit_generation_cap` 必须互斥。本批次 972 条由 EOS 正常结束，28 条达到
生成长度上限；达到上限应保留为警告，不应被静默当作自然完整序列。

## 4. 已完成的基础 QC

仓库已经提供可复用的基础 QC 实现和当前批次结果：

- 配置：`configs/vae_member_a_v1_candidates.yaml`；
- 入口：`experiments/prepare_generated_candidates.py`；
- 实现：`src/gv_eval/generated_candidates.py`；
- QC 表：`data/processed/gv02_03_v2/vae_member_a_v1/candidate_qc.tsv`；
- 通过序列：`data/processed/gv02_03_v2/vae_member_a_v1/qc_pass.fasta`；
- 摘要与清单：同目录下的 `qc_summary.json` 和 `qc_manifest.json`。

当前实现已经检查：

1. FASTA 与生成元数据的 ID、顺序、长度和 SHA-256 一致性；
2. EOS 与长度上限状态互斥；
3. 空序列、非标准氨基酸、50–180 aa 长度范围和完全重复；
4. 残基种类数和单一残基占比定义的低复杂度代理；
5. 与实际训练集及 346 条交接参考的完整序列精确匹配；
6. 相对天然参考平均组成的 Jensen–Shannon divergence，阈值为参考分布 99% 分位数。

当前结果为：1000/1000 通过固定硬规则，174 条有组成偏离警告，28 条有生成上限警告；
没有训练集或参考集精确重合。这里的“QC 通过”仅表示基础序列规则通过，不代表属于 GvpA，
也不代表有功能。

## 5. 复现基础 QC

在仓库根目录运行：

```bash
cd /home/user/wangyuhan/ML-GV02-03
PATH=/home/user/wangyuhan/envs/ml-gv02-03/bin:$PATH \
PYTHONPATH=src \
/home/user/wangyuhan/envs/ml-gv02-03/bin/python -B \
  experiments/prepare_generated_candidates.py qc \
  --config configs/vae_member_a_v1_candidates.yaml
```

命令会先验证冻结输入，再覆盖写入同一批次的可重建 QC 产物。开始扩展 QC 前，建议先运行：

```bash
sha256sum \
  data/generated/sequence_vae/member_a_v1_seed42/vae_candidates.fasta \
  data/generated/sequence_vae/member_a_v1_seed42/generation_metadata.tsv \
  data/generated/sequence_vae/member_a_v1_seed42/generation_manifest.json
```

若哈希与第 2 节不同，应停止并确认是否已切换到新的生成批次，不要把两个批次混合处理。

## 6. 已有家族与相似度证据

本批次已经用成员 A 的竞争性家族分类流水线完成一次暂定分析：

- `candidate_family_classification.tsv`：
  `data/processed/gv02_03_v2/vae_member_a_v1/family/`；
- 逐条联合审计表：
  `data/processed/gv02_03_v2/vae_member_a_v1/candidate_pool_audit.tsv`；
- 当前候选池：147 条 `main_supported_gvpa`、47 条 `ambiguous_exploration`、
  806 条 `excluded`。

联合表已包含 `family_similarity_best_reference`、`family_similarity_effective_identity`、
`family_similarity_query_coverage` 和 `family_similarity_target_coverage`。覆盖率实现已限制在 `[0,1]`。
这些是 BLAST 局部比对的覆盖修正证据，不能冒充全局或半全局比对结果。

如需在家族分类更新后重建三个互斥候选池，运行：

```bash
PATH=/home/user/wangyuhan/envs/ml-gv02-03/bin:$PATH \
PYTHONPATH=src \
/home/user/wangyuhan/envs/ml-gv02-03/bin/python -B \
  experiments/prepare_generated_candidates.py pools \
  --config configs/vae_member_a_v1_candidates.yaml
```

`pools` 必须在 `qc` 和家族分类均完成后执行；不得在缺失家族证据时把 QC 通过者全部放入主池。

## 7. 成员 C 仍需完成的内容

以下内容尚未形成 C 模块的正式交付，不能因本轮已跑出筛选结果而省略：

1. 增加疏水连续片段检查，阈值写入配置，并区分硬失败与警告；
2. 使用广谱结构域库检查 PF00741 以外的额外结构域，而不只是现有 Gvp 家族 profile；
3. 对临界 BLAST 命中执行全局或半全局比对复核，明确输出 `distance_status`；
4. 分别为 147 条当前主池和 47 条歧义池生成距离矩阵、ID 顺序文件和方法说明；
5. 对所有新增阈值记录依据，补充单元测试、输入输出哈希和 Manifest；
6. 若 A 发布新的最终高可信参考，使用新参考完整重跑，不沿用本批次的相似度或组成阈值。

当前已有的
`data/processed/gv02_03_v2/vae_member_a_v1/scoring/candidate_distance.npy`
只服务于现有主分析评分，不是 C 计划中“主池＋歧义池”两套正式距离交付。

建议 C 的新增产物写入独立目录，避免覆盖 B 的冻结输入和当前可复现实验：

```text
data/processed/gv02_03_v2/vae_member_a_v1/c_handoff/
├── candidate_qc_extended.tsv
├── qc_pass_extended.fasta
├── trusted_similarity.tsv
├── main_distance.npy
├── main_distance_ids.json
├── ambiguous_distance.npy
├── ambiguous_distance_ids.json
├── qc_similarity_summary.json
└── qc_similarity_manifest.json
```

建议正式输出至少保留以下稳定字段：

```text
sequence_id
sequence_sha256
qc_pass
qc_reasons
qc_warnings
terminated_by_eos
hit_generation_cap
closest_reference_id
effective_identity
query_coverage
target_coverage
distance_status
```

## 8. 验收标准

C 提交前至少确认：

- 输入恰好 1000 个唯一 ID，且 FASTA、元数据与输出之间没有遗漏或新增；
- 所有输出行均可通过 `sequence_id + sequence_sha256` 回连冻结输入；
- 所有比例字段在 `[0,1]`，距离矩阵为方阵、对称、对角线为 0；
- 距离矩阵 ID 文件与矩阵轴顺序完全一致；
- 主池、歧义池、排除池互斥且并集为 1000；
- 警告与硬失败分开记录，不静默删除候选；
- 新增 QC 和相似度测试通过，同时原有 95 项回归继续通过；
- Manifest 记录配置、参考集、候选、代码和全部正式输出的 SHA-256；
- 报告明确使用的是暂定参考，不把计算证据写成功能确认。

全项目回归命令：

```bash
PYTHONPATH=src \
/home/user/wangyuhan/envs/ml-gv02-03/bin/python -m pytest tests -q -p no:cacheprovider
```

## 9. 交叉审核重点

成员 C 应首先抽查 28 条长度上限警告和 174 条组成偏离警告，再从无警告样本中随机抽取对照。
成员 A 负责审核新增相似度和额外结构域规则是否符合家族背景；成员 D 使用 C 的正式输出时只按稳定
ID 合并，不人工复制结果，也不得用最终筛选排名反向调整 QC 阈值。
