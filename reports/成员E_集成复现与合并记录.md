# 成员 E（庆）：集成、复现与分支合并记录

日期：2026-09-28。对应分工：工程集成、接口核对、端到端测试、复现和主线整理。

## 1. 交付结论

已将 A/B/C/D 的实现和冻结结果整合，修复 Windows 检出与隔离复现问题，新增完整冻结批次复现入口并重写 README。合并后 **304 项测试全部通过，无跳过**；独立 Git 副本中再次通过。

本次对现有 `member_a_v1_seed42` 批次做复现：候选漏斗为 **1000 → 147 → 136 → 114**，三策略 Top-10/20/50 均可重建。参考继续保持暂定状态，名单仍是计算筛选结果。

机器证据位于 [`results/gv02_03_v2/integration_20260928/`](../results/gv02_03_v2/integration_20260928/)；总索引为 `integration_summary.json`，该目录自身的文件哈希在 `manifest.json`。

## 2. 分支与合并范围

| 来源分支 | 合并前提交 | 内容 |
| --- | --- | --- |
| `main` | `236c68a` | A 的参考/家族分类、保守发布包及既有基线 |
| `feat/sequence-vae` | `10d9da8` | B 的 GRU VAE、实际数据划分、1000 条候选和交接 |
| `feat/quality-similarity` | `894f1fa` | C 的 QC、全局相似度、距离矩阵及域/跨膜复核 |
| `feat/v2-evaluation` | `6a54762` | D 的联合门槛、评价/筛选和冻结验收入口 |
| `feat/candidate-quality-control` | `939b3bd` | 已在原主线中的早期 QC，无独立遗漏提交 |

在独立集成分支上逐一合并，保留所有原始提交。冲突发生于 `.gitattributes` 与 README：属性文件同时保留 A 的 LF 规则和 C 的原始字节例外；README 按合并后的实际状态重写。验收完成后将 `main` 快进到集成结果，不重写已有历史；远端成员分支保留用于溯源。

## 3. 实际执行的核对

| 范围 | 本次动作 | 结果 |
| --- | --- | --- |
| A 参考与保守发布 | 核对已有家族产物、输入输出哈希及保守 release | 通过；未重建外部工具模型 |
| B 数据划分与候选 | 核对 train/validation/test、同源簇隔离、词表、候选/元数据/序列哈希 | 277/33/36；1000 条冻结候选一致 |
| C QC 与全局比对 | 独立运行两次；每次 346,000 次参考比对及 10,731/1,081 对池内比对 | 13 个产物两轮字节一致 |
| C 汇总复核 | 重用经哈希验证的 Pfam/TM 原始结果，重新合并 | 915 条复核；3 个发布数据产物等价 |
| D 评价/筛选 | 重算四维评分、九份名单、消融/鲁棒性/随机基线，再独立复跑 | 44 个产物两轮字节一致 |
| 跨环境与 D 发布版对比 | 校验旧哈希、配置、源提交和输入后比较输出 | 42 个数据文件字节一致；2 张 PNG 字节不同 |
| 全仓库测试 | 新建 Python 环境安装依赖后执行 pytest | 304 passed，0 skipped |
| 独立副本 | 独立 Git 克隆；清除 PYTHONPATH，禁用用户包和 pytest 外部插件 | 304 passed，评价复现通过 |
| VAE 冒烟运行 | 真实训练/验证输入、CPU、2 epoch；固定 checkpoint 两次采样，各 32 条 | FASTA、元数据、生成清单字节一致；未评估测试集 |

C 与原发布文件相比，3 个 JSON 仅换行字节不同；其余 10 个产物字节一致。C 汇总的两张 TSV 字节一致，摘要 JSON 内容相同。D 的两张图已检查可读性；数值一致性由原始数据文件证明，不用图片哈希替代科学结果核对。

比较规则在运行前写入实现：同环境要求清单和输出哈希相同；跨环境先验哈希，再核对顺序、ID、整数计数、矩阵形状/NaN 和数值。浮点容差 `rtol=1e-9, atol=1e-12`，变化的 PNG 单列记录。当前 D 的全部数据文件实际上无需浮点容差即可字节匹配。

## 4. 发现并修复的问题

### Windows 旧检出的换行符

首次合并回归有 5 项旧产物测试失败。逐一与 Git blob 和原清单比较，确认旧检出中有 160 个受 LF 规则约束的文本仍为 CRLF；仅在去除 CRLF 后与原提交完全相同的文件上恢复 Git 原始字节。没有改写期望哈希或计算结果。

另有两个历史 TSV 原本就以 CRLF 提交，增加了精确路径的 `-text` 规则，保留其原始字节。C 的审计文件字节例外也完整保留。详见 `checkout_newlines.json`。

### 隔离脚本的平台假设和比较方式

原 `check_frozen_checkout.py` 写死 `venv/bin/python`，Windows 无法执行；且完整 Manifest 比较会把 Python 版本或渲染字节差异混同为科学结果不一致。

现在使用所选解释器，在独立 Git 副本内运行；同步待提交增加/修改/删除，禁用用户包和环境注入，校验原产物哈希，并独立比较科学数据。报告如实说明这是源码隔离，不声称重新安装了依赖。报告写入新路径，保留旧 D 验收文件。

### 单入口与状态说明

新增 `experiments/reproduce_project.py`，串联 A/B 哈希审计、C 两次实算及复核汇总、D 两次评价和发布版比较。它保留固定批次，不覆盖源交接文件；所有新输出写到新目录。

新增 `requirements-integration.txt`，补齐 C 和测试所需的依赖说明。README 列出当前结果、Windows/Linux 命令、模型缺口和成员入口；旧 T05 报告明确保留为历史结果。

## 5. 可执行复现

在 Python 3.11/3.12 虚拟环境中，仓库根目录执行：

```bash
python -m pip install -r requirements-integration.txt
python experiments/reproduce_project.py --output results/reproductions/project_run1
```

全量测试与独立副本检查：

```bash
python -m pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cpu
python -m pytest tests -q
python experiments/check_frozen_checkout.py --report results/reproductions/checkout_acceptance.json
```

本机验收环境：Windows 11、Python 3.12.14、NumPy 2.4.6、pandas 2.3.3、SciPy 1.17.1、Matplotlib 3.11.1、Biopython 1.88、PyTorch 2.8.0+cpu。完整安装记录在 `environment_freeze.txt`。

本机完整重跑产物保存在 `E:/machine/tmp/gv02-integration/project_replay/`。仓库提交紧凑的报告、哈希清单、测试 XML 和冒烟证据；`replay_output_hashes.json` 的路径相对于该重跑目录，不意味着这些重复文件都在验收记录目录中。其他电脑可用上述入口重建。

## 6. 仍需保留的限制

1. **原始权重未交付。** B 的 `best.pt` SHA-256 应为 `d2801a2a130928739f44f9726bcb5bbd7be229fadfe2f05fba5d0f6cc87debdf`。当前仓库没有该文件，不能宣称重放了其 GPU 训练或原始 1000 条采样。此次 CPU 冒烟权重属于独立检查，不替代它。
2. **上游外部工具未全量重算。** A 家族建模、全 Pfam 扫描、PureseqTM 使用原始固定证据，经哈希和 ID 核验；完整重跑需相应工具及模型库。工具说明见 C/A 交接。
3. **科学状态仍暂定。** A 的 14 条注释冲突有 6 条得到注释修正、8 条仍歧义；保守 release 继续使用 346 条。B 训练有低 KL 警告，尚需生成模型改进与独立验证。
4. **没有功能实验证明。** QC、家族支持、排名与人工复核标记分别保留含义。915 条复核中包含 62 条主池候选，工程验收不取消这些警告。

当前完成的是 E 的冻结批次工程验收和主线交付。PPT、答辩稿、真实会议/工时记录及湿实验不属于本次代办范围。
