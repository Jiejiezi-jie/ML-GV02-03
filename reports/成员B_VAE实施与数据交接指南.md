# 成员 B：参考数据交接与序列 VAE 实施指南

## 1. 当前能接收什么，不能把什么当最终版本

成员 A 已提供代码、来源审计、HMM、竞争鉴定、对照和可复现输出；科学验收仍有 14 条 GvpJ 注释冲突。因此当前不能声称“A 已提供无争议最终参考”。

原始 `reference/high_confidence_gvpa.fasta` 有 681 条计算支持序列。为降低冲突传播风险，另建 `data/processed/gv02_03_v2/member_b_handoff/`：

| 文件 | 用途 |
| --- | --- |
| `gvpa_training_eligible.fasta` | 346 条保守筛选参考，只供开发/明确标为暂定的训练 |
| `quarantined_gvpa.fasta` | 335 条暂时隔离参考，不混入训练 |
| `release_audit.tsv` | 每条的来源登录号、SHA-256、簇、保留/隔离理由 |
| `release_manifest.json` | 版本、来源及输出哈希、使用边界 |
| `unresolved_negative_controls.tsv` | 全部 14 条冲突及其原始注释，不是仅保留容易通过的负例 |
| `development_split/` | 新划分的开发用 train/validation/test、词表和 manifest |

隔离规则为：整个混合家族注释簇隔离；14 条冲突所在簇隔离；其最佳 GvpA 比对参考所在簇也隔离。它是保守风险控制，不等于被隔离序列一定错误，也不证明剩下 346 条已实验确认。不得为了凑训练规模把 `possible_gvpa`、负对照、隔离集或生成候选混回训练。

开发可先开始：接口实现、合成样例测试、小批训练、流程联调。正式实验应等 A 发布经核查的新版本，再重新划分、从头训练和生成，不把当前开发 checkpoint 改名成最终模型。

## 2. 数据重划分

运行项目根目录命令：

```bash
python -B experiments/prepare_member_b_handoff.py
python -B experiments/prepare_member_b_split.py --allow-provisional
```

第二个命令默认拒绝非最终参考，`--allow-provisional` 是开发模式声明，不是把状态改成最终。当前脚本先核对 release 的全部输入输出哈希，再从保留下来的参考重新分配训练/验证/测试。沿用 A 已建立的 70% identity、双向覆盖 80% 的完整簇边界，过滤成员后仍保持同簇整体分配；**不沿用 A 的三折编号，也不沿用旧版 1224 条参考的训练划分**。比例目标 80/10/10，种子 42，实际数量按整簇分配记录在 `split_manifest.json`，不能拆大簇强凑比例。

正式版参考改变时，B 需重新聚类或使用新参考全量簇信息；建议同时报告跨集合最邻近 identity/覆盖分布。CD-HIT 是启发式聚类，“簇 ID 不重叠”不保证任意跨集合序列都不相似。禁止把 VAE 测试集用于早停、参数选择、生成温度选择或筛选 checkpoint。

A 的家族 profile 可能使用了与 B 验证/测试序列同源的审校锚点。家族通过率应作为固定注释工具的代理指标，不能叫“独立分类准确率”。B 的生成模型训练与测试仍需自身严格隔离。

当前实际划分为 31 簇、训练 254 条、验证 60 条、测试 32 条（合计 346）。集合较小且有大簇，建议优先做小模型和重复种子实验；切勿通过拆簇人为修正比例。验证/测试指标应附样本量，不作过强泛化结论。

## 3. 现有代码的具体问题

`gv/design/src/` 是参考实现，当前不能直接作为正式训练器。检查发现：

| 位置 | 现状 | B 必须处理 |
| --- | --- | --- |
| `data_loader.py` 词表 | AA 0–19、gap 20、PAD 21、SOS 22、EOS 23 | 与本项目 PAD 0、BOS 1、EOS 2、UNK 3、AA 4–23 不兼容；全链路只加载同一份 JSON |
| `vae_gvp.py` decode/generate | 使用 21 作起始 token | 这是旧词表 PAD，不是 SOS；禁止硬编码 |
| `vae_loss()` | `x != 23` 作 mask | 屏蔽了 EOS，未屏蔽 PAD；改为忽略 PAD，EOS 应计入重构损失 |
| `generate()` | 任一条 EOS 即停止全 batch | 用逐样本 `finished`，直到全部完成或达到上限 |
| `generate()` | 达目标长度后人为补 EOS | 不得把补写 EOS 记成模型自然结束；保存 cap/stop_reason |
| `_clean_sequence()` | 删除特殊 token、截断重复片段 | 会抹除原始生成问题；保存 raw tokens 后单独 QC，不静默修饰 |
| 编码器 | 对含 PAD 的整段 LSTM 取末状态 | 使用真实长度和 packed sequence 或正确 mask，避免 PAD 影响潜变量 |
| 长度头 | `.mean().item()` 决定 batch 循环长度 | 截断梯度且 batch 样本互相影响；基础版先用固定最大步数+EOS，若保留长度头需直接监督其输出 |
| 长度损失 | 只依赖目标真实长度 | 对长度预测头无有效训练信号，不应称为已训练的长度控制 |
| `load_state_dict()` | 自动补缺失参数且 `strict=False` | 正式恢复必须检查架构/词表/键名，缺失即报错 |
| 数据加载 | 从文件名或 header 推断链/物种 | `train.fasta` 和 `nat_<hash>` 不含可信物种；hash 后缀会被误作物种，造成伪条件或空集 |
| `train.py` | 仅训练循环和周期性 state_dict | 缺验证、早停、best checkpoint、曲线、完整恢复状态和元数据 |

推荐第一版使用 `configs/generation.yaml` 指定的无条件 GRU 序列 VAE；若选现有条件 LSTM，必须明确变更配置和报告。当前 346 条数据不适合为每个序列 ID 建一个物种条件，先用无条件模型或统一 unknown。不要把目录路径交给会同时读取 train/validation/test 的加载器；显式传入三个文件。

## 4. VAE 实施基线

建议配置：embedding 32、hidden 128、latent 32、Adam lr 0.001、batch size 64、最多 100 epochs、KL beta 0.1、warm-up 20 epochs、early stopping patience 15。这些是起始建议，不是已训练最优超参数。

编码：`BOS + amino_acids + EOS`，再 PAD 到 batch 最大长度。AA 上限 180 与 token 上限 182 分开命名；长度恰好为 180 的真实序列不能被 `>= max_len-2` 静默丢弃。遇到非法字母、空序列、超长序列应明确报错或记入剔除表。

训练解码：输入 `x[:, :-1]`，目标 `x[:, 1:]`。序列结束 EOS 仍是有效目标；只有 PAD 屏蔽。重构损失按有效 token 总数归一，KL 按每条序列归一并说明 beta 的尺度。单独记录 recon、KL、beta、total，监测 posterior collapse；不能只画 total。

验证：`model.eval()` 和 `no_grad()`；固定验证采样策略，推荐用 `z=mu` 作确定性重构，另记采样 ELBO。早停指标用固定 beta 的 validation loss 或固定重构 NLL，不能把 warm-up 前后不同 beta 的 total loss 混作同一标尺。验证按 token 总量累积再除，避免长短 batch 权重错误。

保存 `best.pt` 与 `last.pt`：包含 model/optimizer/scheduler 状态、epoch、best metric、patience counter、Python/NumPy/Torch CPU/CUDA RNG 状态、架构、完整配置、词表及哈希、release/split 哈希、软件版本、设备。恢复后校验词表、维度和输入数据哈希，不能静默适配旧权重。

## 5. 曲线和模型评价

每 epoch 写 `training_history.tsv`：epoch、train/val reconstruction NLL、KL、固定 beta validation loss、当前 beta、learning_rate、非 PAD token accuracy、epoch_seconds、best_epoch 标志。导出训练/验证曲线和 KL 曲线；不要仅截取终端输出。

模型选定后只进行一次最终 test 评价，报告 token NLL/accuracy、精确序列重构比例、长度分布、KL 和潜变量方差。低 loss 不代表生成有功能。再报告生成有效率、自然 EOS 比例、cap 比例、唯一率、精确训练重合率、氨基酸组成、长度分布、低复杂度比例、统一家族鉴定通过率。去重后统计与全部 1000 条统计都要给出，不能只展示筛选后的好样本。

## 6. 生成 1000 条与逐序列元数据

使用锁定的 best checkpoint、固定 prior `z~N(0,I)`、种子 42、temperature 1.0（默认 top_k=0、top_p=1.0）。输出 1000 个原始采样尝试，ID 如 `vae_seed42_000001`。若最终合格不足 1000，报告不足，不静默补采来隐藏失败；额外采样批次另立记录。

每条保存以下字段，TSV 一行对应 FASTA 一条，原始 token 另存 JSONL：

| 字段 | 约定 |
| --- | --- |
| sequence_id / sample_index / batch_id | 唯一稳定；sample_index 明确从 0 或 1 开始 |
| generator_type | `sequence_vae` |
| model_checkpoint / checkpoint_sha256 | 实际使用的 best checkpoint |
| reference_release_id / split_manifest_sha256 / vocabulary_sha256 | 数据与词表可追溯 |
| generation_seed / sample_seed / latent_id | 可复现随机流及潜变量索引；可保存 latent 文件和哈希 |
| temperature / top_k / top_p | 实际采样参数 |
| max_amino_acids / raw_token_length / sequence_length | 分别记录残基上限、原始 token 数、输出残基数 |
| terminated_by_eos | 仅模型实际采到 EOS 才为 true |
| hit_generation_cap / stop_reason | `eos`、`length_cap`、`invalid_token` 等互斥原因 |
| sequence_sha256 | 解码后大写 AA 序列 SHA-256 |
| device / torch_version / sampling_batch_size | 影响随机性和数值行为的运行信息 |

EOS 规则：每条第一个 EOS 后的内容不再计入序列；完成样本不影响其他样本。生成时屏蔽 PAD/BOS/UNK 的采样，只允许 20 种 AA 和 EOS。若恰好在最后允许步采到 EOS，按自然终止记录；未采到 EOS 才记 cap。人为补 EOS 不改变 cap 状态。空序列也保留 raw attempt 和原因，不伪造残基填充。

本项目 QC 的 `terminated_by_eos`、`hit_generation_cap`、`generation_cap_warning` 等含义必须与元数据一致。默认 180 AA 上限时，需要允许在 180 AA 后再作 EOS 判断，明确算法顺序并测试边界。

## 7. 随机种子与复现验收

设置 Python random、NumPy、Torch CPU/CUDA、DataLoader generator/worker seed，关闭 cuDNN benchmark，按需启用 deterministic algorithms；固定依赖、设备、batch size、采样顺序。GPU 若无法逐字节确定，明确报告范围。

必须测试：同 checkpoint+同配置+同设备+同种子，两次生成 1000 条的 FASTA 与元数据关键列一致；不同种子产生差异；checkpoint 保存/恢复前后 logits 一致。若使用全局随机流，batch size 改变可能改变结果；不要宣称跨 batch size 复现。逐样本 seed 或预先冻结潜变量/采样流可降低这种影响。

关键单元测试：PAD 不计损失、EOS 计损失、teacher-forcing 移位正确、不同结束时刻的 batch、180 AA 边界、无 EOS cap、空输出、非法 token、严格 checkpoint、train/val/test ID+序列哈希+簇互斥。另加 8–16 条序列的小集过拟合测试，以检查模型确实能学会重构。

## 8. 目录与交付清单

建议 B 新建 `src/gv_eval/generation.py`、`experiments/train_generator.py`、`experiments/generate_candidates.py` 和对应 tests。模型置 `models/generator/sequence_vae/<release>/<run>/`；1000 条新候选及元数据置 `data/generated/sequence_vae/<batch>/`；曲线和评价置 `results/gv02_03_v2/generation/<run>/`。这些是建议新增的接口，当前没有假装它们已实现或已有权重。

正式交付检查：最终参考 release 已明确；新簇划分和词表固定；训练/验证/早停记录齐全；best/last checkpoint 可严格恢复；曲线和一次性 test 评价齐全；1000 条原始候选及 metadata/raw tokens 全覆盖；EOS/cap 边界测试通过；固定种子复现有文件证据；生成结果进入 A 家族判定及统一 QC，而不是直接宣称为 GvpA。

当前家族入口 `run_family_classification.py --candidates <新候选>` 会重建参考并覆盖配置指定目录。对新批次先复制 `configs/family_classification.yaml`，将 `reference_dir`、`profile_dir`、`result_dir` 改成该批次独立目录，避免覆盖 A 冻结证据；候选 TSV 写在该 `reference_dir` 的父目录，目录设计也要避免重名。

GAN 是可选项，只有 VAE 主流程通过验收且时间允许再做。GAN 必须使用同一最终数据划分、相同原始采样预算、相同 QC 和家族工具，报告模式坍塌/重复率。不能为 GAN 单独放宽门槛以制造优势。
