# 成员 B 的 VAE 最小参考包

整理日期：2026-09-22。本目录只保留补充材料中与 B 的 VAE 开发直接相关的 6 个文件。
未包含其他模型、整包数据、历史图表、评分器、论文、识别模块或基因线路模块。

**已有真实 VAE 源码和训练入口，没有训练好的权重；仍需修正后接入本项目 V2。**

## 1. 保留文件

| 文件 | 作用 |
| --- | --- |
| 本说明 `README.md` | B 的入口、已知问题与接入要求 |
| [design/src/models/vae_gvp.py](design/src/models/vae_gvp.py) | `GVAE` 模型、`vae_loss` 和 `GVAE.generate()` 生成方法 |
| [design/src/models/__init__.py](design/src/models/__init__.py) | 模型模块包入口 |
| [design/src/data_loader.py](design/src/data_loader.py) | `GvpDataset`、旧词表、批处理和 `decode()` |
| [design/src/train.py](design/src/train.py) | 已精简为只训练 VAE，保留原 `train_vae()` 算法 |
| [design/requirements-vae.txt](design/requirements-vae.txt) | PyTorch、NumPy、Biopython 三项直接依赖 |

训练脚本已去掉其他模型的导入和训练分支，且要求显式提供 `--fasta`，避免默认混读多家族目录。
模型和数据加载核心源码未做算法修改。生成入口在 `GVAE.generate()` 内；
B 后续应封装独立生成程序，结合 `GvpDataset.decode()` 导出 FASTA 和生成元数据。

## 2. 使用现有数据和配置

- [本项目 V2 配置](../configs/generation.yaml)：计划中的 GRU-VAE 参数和候选数量。
- [暂定训练/验证/测试划分](../data/processed/gv02_03_v2/generator_data/)：可用于开发准备，
  正式训练仍需等待 A 的高可信 GvpA 参考并重新划分。
- [已有 GvpA 原始来源](../data/raw/design/GvpA.fasta)：补充目录原有的同名 856 条序列与该文件
  规范化换行后内容相同，因此未重复保留。

在安装好三项依赖的环境中，从仓库根目录可用以下命令查看训练参数：

```bash
python gv/design/src/train.py --help
```

这不是正式训练命令。训练时必须显式传入核验后的 GvpA 数据；旧数据加载器会按文件名和标题
推断家族及物种，接入 `train.fasta` 等规范化 ID 数据时需要调整，避免过滤为空或把 ID 当作物种。
权重将写入运行工作目录的 `checkpoints/vae_epoch*.pt` 和 `checkpoints/vae_final.pt`。

## 3. 与现有 V2 方案的差异

| 项目 | 这份参考代码 | 现有 V2 计划 |
| --- | --- | --- |
| 网络 | 条件 LSTM-VAE | GRU-VAE |
| embedding / hidden / latent 默认维度 | 128 / 256 / 64 | 32 / 128 / 32 |
| 特殊符号 | PAD=21、SOS=22、EOS=23，另有 gap | PAD=0、BOS=1、EOS=2、UNK=3 |
| 验证与训练控制 | 尚无验证集、早停、KL warm-up | 配置已规划早停和 KL warm-up |
| 生成元数据 | 尚不完整 | 需保存种子、权重哈希、EOS、截断等 |

B 应统一模型结构、词表和配置；不能直接混用两套编号和权重。
旧训练入口不会自动读取 `configs/generation.yaml`。

## 4. 训练前需要修正的问题

以下为静态检查发现的问题，本次仅精简材料，没有修复模型算法：

1. `data_loader.py` 定义 SOS=22，但模型 `decode()` 和 `generate()` 用 PAD=21 作为起始输入。
2. `vae_loss()` 屏蔽 EOS=23，没有屏蔽真正的 PAD=21；需核对重构目标和填充位置。
3. 生成时任意样本出现 EOS 就停止整个 batch；需逐样本区分正常结束、截断和清洗。
4. 长度预测经 `.item()` 转整数控制循环，不能经此路径获得梯度；长度损失又只由真实目标计算。
5. 兼容加载会自动补齐缺失参数并使用 `strict=False`，正式恢复时需验证权重完整性。

还需完成固定种子、按同源簇划分的数据接入、验证/早停、KL warm-up、训练曲线和模型状态记录。
生成结果应交由本项目统一 QC 与评价流程处理，不能直接视作已确认的功能 GvpA。

## 5. 本次检查边界

已核对精简包的 Python 语法、本地导入依赖和文档链接；保留的模型、数据读取源码与导入版本一致，
训练函数逻辑未变。依赖列表尚未经过运行环境版本锁定。
当前检查环境没有 PyTorch，因此未验证前向/反向计算、正式训练或候选生成。
