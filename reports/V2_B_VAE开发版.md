# V2 成员 B：序列 VAE 开发版

## 完成范围

本次将组员提供的条件 LSTM-VAE 最小参考包接入项目工程，形成与
`configs/generation.yaml` 一致的无条件 GRU 序列 VAE 开发版。当前只完成代码、测试和
暂定数据冒烟验证，没有将暂定数据训练结果作为正式模型，也没有生成正式1000条候选。

实现内容：

- 统一词表为 `PAD=0、BOS=1、EOS=2、UNK=3` 加20种标准氨基酸；
- 使用 padding-aware 双向 GRU 编码器、重参数化潜变量和 GRU 自回归解码器；
- 重构损失只忽略 PAD，EOS 作为正常预测目标；
- 使用 KL warm-up、梯度裁剪、验证集早停和 posterior-collapse 警告；
- checkpoint 保存模型、优化器、随机状态、配置和输入哈希，支持严格恢复；
- 每条生成序列独立处理 EOS，不因同批其他序列结束而提前截断；
- 生成时禁止 PAD、BOS 和 UNK，并记录 checkpoint 哈希、种子、采样参数、潜变量编号、
  EOS、长度上限、长度和序列哈希；
- 暂定数据默认禁止训练和生成，开发运行必须显式使用 `--allow-provisional-data`。

## 开发命令

下面的命令只能用于当前暂定数据上的代码验证：

```bash
python -B experiments/train_generator.py \
  --config configs/generation.yaml \
  --device cpu \
  --allow-provisional-data \
  --maximum-epochs 1 \
  --output-dir /tmp/gv02_vae_smoke/model

python -B experiments/generate_candidates.py \
  --config configs/generation.yaml \
  --device cpu \
  --allow-provisional-data \
  --candidate-count 5 \
  --checkpoint /tmp/gv02_vae_smoke/model/best.pt \
  --output-dir /tmp/gv02_vae_smoke/generated
```

正式运行不应使用 `--allow-provisional-data`。成员 A 完成高可信 GvpA/GvpJ 竞争性核验后，
应更新配置中的数据状态、按同源簇重新划分，并从头正式训练。

## 验证结果

- 自动化测试：64项全部通过，其中VAE新增测试6项；
- 暂定真实数据：训练978条、验证123条、测试123条；
- CPU冒烟训练：1个epoch成功完成，约9秒；
- checkpoint、训练历史、曲线和Manifest均成功生成；
- 从冒烟checkpoint生成5条候选，5条均独立产生EOS，没有命中180 aa长度上限；
- 上述数值只验证工程链路，不用于评价模型质量。

## 未完成和依赖

正式实验仍依赖成员 A 提供：

```text
high_confidence_gvpa.fasta
reference_classification.tsv
GvpA/GvpJ竞争性核验结果
```

收到后需重新划分数据、在正式CUDA/PyTorch环境中训练至早停、检查验证KL及生成质量，
随后生成1000条候选及 `generation_metadata.tsv` 并交给成员 C 做完整QC。
