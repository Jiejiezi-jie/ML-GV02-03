# GV02-03：GvpA 候选序列的生成、评价与筛选

课程项目 Track 1 / GV02-03。项目使用序列 VAE 生成蛋白候选，结合质量控制、GvpA/GvpJ 竞争性家族鉴定和四维评价，比较等权加权、Pareto、分维度轮转三种筛选策略。

**当前主线：`member_a_v1_seed42` 暂定参考批次。** A/B/C/D 的实现及交接已汇集到 `main`；E 完成合并后测试与冻结批次复现。这里的名单表示计算筛选结果，尚无结构、表达、组装或功能实验确认。

## 1. 已完成的结果

| 筛选步骤 | 数量 | 说明 |
| --- | ---: | --- |
| A 交给 B 的保守参考 | 346 | 仍为 `provisional_not_scientifically_final` |
| 训练 / 验证 / 测试 | 277 / 33 / 36 | 按 31 个同源簇隔离；对应 16 / 9 / 6 簇 |
| B 冻结 VAE 候选 | 1,000 | 全部唯一；972 条自然 EOS 结束，28 条触及生成上限 |
| 基础 QC 通过 | 1,000 | 通过格式、长度等规则不等于家族或功能确认 |
| GvpA 家族支持 | 147 | 另有 47 条 A/J 歧义、806 条排除 |
| QC + 家族 + 域门槛通过 | 136 | PF00741 分数 ≥25 bits、模型覆盖 ≥0.95 |
| 四维证据足够、可排名 | 114 | 其余 22 条评分证据不足，保留在审计表 |
| 三种策略 | 9 份名单 | 每种输出 Top-10、Top-20、Top-50 的 TSV 和 FASTA |

已完成逐维消融、200 次权重扰动、27 组域/配对覆盖门槛组合、1,000 次同池随机基线，以及与 B 原局部比对版本的名单对比。

数据依据：[D 结果摘要](results/gv02_03_v2/d_evaluation_member_a_v1/summary.json)、[训练记录](models/generator/sequence_vae/member_a_v1_seed42/training_manifest.json)、[生成记录](data/generated/sequence_vae/member_a_v1_seed42/generation_manifest.json)。

## 2. 主流程与评价规则

```text
A：参考来源审计、GvpA/GvpJ 竞争核验、保守参考发布
  → B：按同源簇划分 → GRU 序列 VAE → 冻结 1000 条候选和生成元数据
  → C：质量控制 → 全局参考比对 → 同类距离矩阵 → 域/跨膜预警
  → D：QC、家族、域联合门槛 → 四维评分 → 三策略 Top-K → 鲁棒性分析
  → E：接口/哈希校验 → 两次运行比较 → 独立副本测试 → main 集成
```

四维含义：

- **约束满足度**：PF00741 域模型的连续支持分数；该域同时覆盖部分 GvpJ，不能单独证明 GvpA 身份。
- **保守性**：相对于暂定 GvpA 参考的 23 个统计保守位点的一致程度。
- **新颖性**：与可靠天然参考的最近全局距离。
- **候选独特性 / 集合多样性**：在同类、联合门槛合格池中比较候选距离；集合另外报告配对覆盖率和未知距离的上下界。

未解析距离保留为 NaN，不能补为 1。排名要求可靠配对比例至少 50%，并检查 25%/50%/75% 的敏感性。人工复核警告保留在名单中；不从歧义或排除池补足 Top-K。

## 3. 快速复现评价与筛选

需要 Git 和 **Python 3.11 或 3.12**。使用完整克隆，保留历史提交：冻结入口按 B 的 `10d9da8` 与 C 的 `894f1fa` 读取 Git 原始字节。ZIP 下载、仅复制源码或缺失历史对象的浅克隆不能直接使用该入口。

```bash
git clone https://github.com/Jiejiezi-jie/ML-GV02-03.git
cd ML-GV02-03
```

Windows PowerShell：

```powershell
py -3.12 -m venv .venv
.venv/Scripts/python.exe -m pip install -r requirements-frozen-v2.txt
.venv/Scripts/python.exe experiments/run_full_experiment.py --output results/reproductions/evaluation_run1
```

Linux / WSL：

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-frozen-v2.txt
.venv/bin/python experiments/run_full_experiment.py --output results/reproductions/evaluation_run1
```

该命令核验冻结输入，重新评分、筛选和运行实验，再独立复跑比较全部 44 个产物。**输出目录必须是新的或空目录**；第二次手动运行请换目录。此模式无需 GPU、Torch 或原始 checkpoint。

`--verify-existing` 用于相同代码和环境下严格比较已有运行，包含版本与图片哈希；跨环境核对发布结果请用下一节的 E 入口。

## 4. E 的完整冻结批次复现

在上一节的环境中安装扩展依赖，再执行（以下 `python` 指虚拟环境中的解释器）：

```bash
python -m pip install -r requirements-integration.txt
python experiments/reproduce_project.py --output results/reproductions/project_run1
```

此入口会：

1. 核验 A 参考/家族产物、保守发布包、B 的实际划分与冻结输入哈希。
2. 两次重算 C 的基础/扩展 QC、346,000 次候选—参考全局比对及主池/歧义池距离。
3. 结合已发布且经哈希核验的 Pfam / PureseqTM 表，重新汇总 C 复核名单。
4. 两次运行 D 的四维评分、三策略筛选和鲁棒性实验。
5. 与成员发布产物核对，保存 `reproduction_report.json` 和全目录 `manifest.json`。

D 继续读取冻结的 B/C 提交；C 的本轮重算用于独立核对其等价性。改动当前工作区的数据不会自动创建新的研究批次。

同一环境的两轮产物要求字节一致。跨环境先校验原始文件哈希，再核对 ID、顺序、整数计数、NaN 位置及科学数值；浮点容差为 `rtol=1e-9, atol=1e-12`，PNG 字节差异单列记录。

**2026-09-28 实跑：** C 的 13 个产物两轮字节一致；D 的 44 个产物两轮字节一致。与 D 发布结果相比，42 个数据文件字节一致，仅 2 张 PNG 渲染字节不同。详见 [E 复现记录](reports/成员E_集成复现与合并记录.md)及[机器可读证据](results/gv02_03_v2/integration_20260928/)。

## 5. 测试与模型重放

全仓库测试还需要 Torch；本次使用 CPU 版 PyTorch 2.8.0：

```bash
python -m pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cpu
python -m pytest tests -q
python experiments/check_frozen_checkout.py --report results/reproductions/checkout_acceptance.json
```

本次 **304 项测试通过，无跳过**，独立 Git 副本中再次通过。隔离检查使用当前解释器及已安装依赖，清除 `PYTHONPATH`、禁用用户包和 pytest 外部插件，并校验科学产物；它不是另一次依赖安装。

### 原始 checkpoint 的缺口

B 的 `best.pt` 没有提交到 Git。原记录的 SHA-256：

```text
d2801a2a130928739f44f9726bcb5bbd7be229fadfe2f05fba5d0f6cc87debdf
```

因此上述验收不能证明 B 原训练或原始 1000 条采样已被重放。需从 B 取得匹配权重，并核对其环境：Python 3.10.20、Torch 2.8.0+cu128、CUDA 12.8、NumPy 1.26.4；跨设备采样不保证逐字节一致。

为了检查模型入口，E 已用真实训练/验证划分完成 **2 epoch CPU 冒烟训练和两次 32 条采样**；两次候选、元数据和生成清单字节相同，未评估测试集，也未替换正式冻结输入。重跑该检查：

```bash
python experiments/train_generator.py --config configs/generation_member_a_v1.yaml --device cpu --allow-provisional-data --maximum-epochs 2 --skip-test-evaluation --output-dir results/reproductions/vae_smoke
python experiments/generate_candidates.py --config configs/generation_member_a_v1.yaml --checkpoint results/reproductions/vae_smoke/best.pt --device cpu --allow-provisional-data --candidate-count 32 --output-dir results/reproductions/smoke_candidates
```

A 的外部家族建模、C 的全 Pfam 扫描与 PureseqTM 本轮只核验已有证据。独立重建这些步骤所需的工具/模型说明见下方交接文档。

## 6. 从哪里读结果

主要结果位于 [`results/gv02_03_v2/d_evaluation_member_a_v1/`](results/gv02_03_v2/d_evaluation_member_a_v1/)：

| 文件 | 用途 |
| --- | --- |
| `candidate_audit.tsv` | 全部 1000 条的家族、QC、域、排名资格与排除原因 |
| `eligible_pool_scores.tsv` | 136 条联合门槛合格候选；保留评分不足原因 |
| `ranking_pool_scores.tsv` / `ranking_pool.fasta` | 114 条可参与四维排名的候选 |
| `{weighted_sum,pareto,dimension_round_robin}_top{10,20,50}.{tsv,fasta}` | 三种策略的有序名单及序列 |
| `ranking_distance.npy` / `ranking_distance_ids.json` | 按 ID 对齐的距离矩阵与缺失距离 |
| `strategy_summary.tsv` / `strategy_overlap.tsv` | 策略指标、复核数量、多样性与名单重叠 |
| `ablation_*` / `weight_robustness.tsv` / `threshold_coverage_sensitivity.tsv` / `random_*` | 消融、权重、阈值与随机基线实验 |
| `summary.json` / `manifest.json` | 摘要、源提交、配置和输入输出 SHA-256 |

## 7. 分工、分支与交接

| 角色 | 来源 | 入口说明 |
| --- | --- | --- |
| A：参考与家族核验 | 原 `main`，`236c68a` | [A 工作过程](reports/成员A_完整工作过程与交接记录.md)；[方法与验收](reports/成员A_家族鉴定方法与验收.md) |
| B：VAE 与候选 | `feat/sequence-vae`，`10d9da8` | [B→C 交接](reports/成员B_to_C_候选序列交接说明.md)；[训练实验](reports/V2_B_VAE完整实验.md) |
| C：QC、相似度与预警 | `feat/quality-similarity`，`894f1fa` | [C 统一交接](docs/c-final-handoff.md)；[Pfam](results/gv02_03_c_domains_member_a_v1/README.md)；[跨膜](results/gv02_03_c_tm_member_a_v1/README.md) |
| D：四维评价与筛选 | `feat/v2-evaluation`，`6a54762` | [D 交接](docs/D_现有批次评价与筛选交接.md) |
| E（庆）：集成、复现与整理 | 本次 `main` 集成 | [E 验收记录](reports/成员E_集成复现与合并记录.md) |

`feat/candidate-quality-control` 的 `939b3bd` 早已包含在主线。各功能分支的提交历史全部保留；原远端分支保留作协作记录，使用者直接从 `main` 获取完整工程。

主要目录：`src/gv_eval/` 是实现，`experiments/` 是命令入口，`configs/` 是参数，`data/generated/` 保存 B 的冻结输入，`data/processed/gv02_03_v2/` 保存参考和候选证据，`results/` 保存结果，`reports/` 和 `docs/` 保存交接与方法。

## 8. 结论边界与历史材料

- 参考仍是暂定发布。A 的 14 条注释冲突中，6 条经数据库复核为旧标签问题，8 条仍歧义；没有据此自动恢复隔离簇。见[复核记录](results/gv02_03_v2/family/adjudication/summary.json)。
- B 的训练记录存在潜变量 KL 偏低警告（posterior collapse）；生成唯一序列和通过 QC 均不能证明生成模型已学到有效功能空间。
- C 的 915 条人工复核标记可与家族/评分池重叠；主池中 62 条有复核标记。预警不会被“测试通过”消除。
- 旧 200 条 T05 候选及 `results/gv02_03/`、`results/gv02_03_qc/`、M1–M5 报告和 PPT 保留作历史基线。旧报告里的数字不代表当前 VAE 批次。
- 旧流程入口为 `python experiments/run_full_experiment.py --workflow legacy --config configs/gv02_03.yaml`；需 Linux/WSL 的 HMMER、MAFFT、BLAST、CD-HIT 等工具，参考 `environment.yml`。本次没有重新验收该完整外部工具环境。
- `gv/` 为早期 VAE 参考代码；当前模型实现是 `src/gv_eval/vae.py` 与 `generation.py`。

当前已完成冻结批次的工程集成与复现。后续研究应补交原权重、解决参考身份歧义，并在新增正式参考或候选批次上完整重跑；功能结论仍需实验。
