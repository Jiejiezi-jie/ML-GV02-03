# C 组非目标结构域与域架构复核（member A v1）

这是对冻结 B→C 批次的独立、预警性检查，不改写 B 的基础 QC、家族分类或候选池，也不把模型命中解释为实验验证。全部序列和计算均在本机。

## 方法

- 输入：1,000 条 VAE 候选和 346 条暂定 GvpA 参考；文件 SHA-256 见 `manifest.json`。
- 用 PyHMMER 0.12.3 将所有序列与 [EMBL-EBI Pfam-A](https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/) 的 30,134 个 profile HMM 比对，使用每个模型自带的 gathering（GA）序列及域阈值。模型压缩包 MD5 为 `7ab3c4e215d0daaea3004e37c4e24f8a`。
- `PF00741` 是目标 Gas_vesicle 模型。其他显著 Pfam 命中分为与 PF00741 不重叠的非目标命中，以及存在任意重叠的替代家族命中；两者只触发人工复核，不自动认定为“新增功能域”。
- 若一条候选的 **PF00741 命中片段数**超过参考序列中观察到的最大值，也触发人工复核。多片段可能来自类重复序列，也可能是 HMMER 将一段复杂序列拆成相邻或重叠命中；**不能据此直接宣称存在两份完整结构域**。逐条证据保留模型坐标和覆盖比例。

## 结果

346 条参考序列均只检出 1 个 PF00741 命中片段，没有其他显著 Pfam 命中。1,000 条候选中同样**没有检出非 PF00741 的显著 Pfam 命中**。

| 原候选池 | 总数 | 单段 PF00741、无额外域证据 | 无显著 Pfam 命中 | 多段 PF00741，需复核 |
| --- | ---: | ---: | ---: | ---: |
| main_supported_gvpa | 147 | 147 | 0 | 0 |
| ambiguous_exploration | 47 | 46 | 0 | 1 |
| excluded | 806 | 382 | 395 | 29 |

共 30 条多段 PF00741 候选需要人工复核，其中模糊池的 `vae_candidate_000113` 为 13–51 和 53–86 位；其余 29 条均在原排除池。完整名单、所有候选的无命中状态和逐段模型覆盖证据分别见 `candidate_domains.tsv`、`candidate_domain_hits.tsv`。参考对照在对应的 `reference_*.tsv`；`summary.json` 与 `manifest.json` 给出统计和哈希。

“无显著 Pfam 命中”不等于证明没有结构域；Pfam 也不覆盖所有结构、重复和膜相关特征。本检查不能替代 B 已做的竞争性家族鉴定、独立结构预测或实验验证。此前的跨膜预测是[另一项独立预警检查](../gv02_03_c_tm_member_a_v1/README.md)。

## 本地复现

在仓库根目录，准备同一冻结 B 快照以及已校验的 Pfam-A 压缩模型；先在本机环境安装 `pyhmmer==0.12.3` 与 `psutil==7.2.2`，再运行：

```powershell
$env:PYTHONPATH = (Resolve-Path src).Path
.venv/Scripts/python.exe -m gv_eval.c_domains `
  --config configs/c_handoff_member_a_v1.json `
  --handoff-root .venv/b-handoff-10d9da8 `
  --pfam .venv/c-tools/pfam/Pfam-A.hmm.gz `
  --output results/gv02_03_c_domains_member_a_v1_rerun `
  --cpus 4
```

脚本直接流式读取 `.gz`，无需解压或预处理模型库；输出目录必须为空，不覆盖现有结果。
