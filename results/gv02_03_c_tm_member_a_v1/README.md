# C 组跨膜区本地复核（member A v1）

本次只对冻结的 B→C 交接数据增加一项**预警性**拓扑检查；不修改 B 的 QC、家族判断和候选池，也不据此宣称实验验证。所有序列留在本机，没有上传到线上服务或提供邮箱。

## 方法与判定

- 输入：1,000 条 VAE 候选，以及 346 条 B 交接的 GvpA 合格参考序列。输入 SHA-256 见 `manifest.json`。
- 工具：[PureseqTM](https://github.com/PureseqTM/PureseqTM_Package)，固定源码版本 `4bb463f5b973f6b7b48bdab55b5cfd150f1741e3`，Linux/WSL 本地运行，快速模式 `-m 0`。逐条核验预测标签只含 `0/1`，且长度与原序列相同。
- 把连续的 `1` 解释为一个预测跨膜片段；区间为从 1 开始、两端都包含的氨基酸位置。346 条参考序列的预测片段数全为 0。因此候选只要出现至少 1 段，就标为 `extra_tm_manual_review`；这是相对于这批参考的异常信号，不是自动淘汰条件。

## 结果

| 候选池 | 总数 | 新增跨膜复核警告 |
| --- | ---: | ---: |
| main_supported_gvpa | 147 | 0 |
| ambiguous_exploration | 47 | 0 |
| excluded | 806 | 7 |

7 条需人工复核的候选及预测区间：

| 序列 ID | 区间（aa） | 片段数 |
| --- | --- | ---: |
| vae_candidate_000067 | 49–67；80–92 | 2 |
| vae_candidate_000101 | 21–39 | 1 |
| vae_candidate_000365 | 10–28 | 1 |
| vae_candidate_000464 | 24–43 | 1 |
| vae_candidate_000600 | 17–36 | 1 |
| vae_candidate_000613 | 21–41 | 1 |
| vae_candidate_000636 | 37–56 | 1 |

逐条结果在 `candidate_tm.tsv`，参考分布在 `reference_tm.tsv`，机器摘要和输入/工具/输出哈希分别在 `summary.json`、`manifest.json`。模型的 2 状态预测并不等于真实膜定位；尤其短片段和边界需要结合更多证据人工复核。“无新增警告”也不等于确认没有异常膜结构。

## 复现

在带有冻结 B 快照的仓库根目录，使用 Ubuntu/WSL、Python 3 和上述版本的 PureseqTM 源码（确保脚本与模型文件为 LF 换行）运行：

```bash
PYTHONPATH=src python3 -m gv_eval.c_tm \
  --config configs/c_handoff_member_a_v1.json \
  --handoff-root .venv/b-handoff-10d9da8 \
  --tool .venv/c-tools/PureseqTM-linux \
  --cache .venv/c-tools/tm-run-v1 \
  --output results/gv02_03_c_tm_member_a_v1_rerun \
  --workers 4
```

`--cache` 只保存本地逐序列预测中间件；脚本使用已有结果时仍会重新核验每条预测的序列和标签格式。请指定空的新输出目录，脚本不会覆盖已有结果。此结果与先前 C 主报告相互独立，未重新排序或更改候选池。
