
# 成员 A：参考集与竞争性家族鉴定

## 运行与交付入口

在项目根目录、具备 Python、PyYAML、HMMER、MAFFT、CD-HIT、BLAST+ 的环境执行：

```bash
python -B experiments/run_family_classification.py                      # 1 分类与参考产物
python -B experiments/verify_family_deliverables.py                     # 2 生成验收表与待核查对照表
python -B experiments/adjudicate_family_conflicts.py                    # 3 按 NCBI 快照裁决冲突
python -B experiments/verify_family_deliverables.py                     # 4 用裁决结果重算验收
python -B experiments/prepare_member_b_handoff.py                       # 5 冻结保守发布包
python -B experiments/prepare_member_b_split.py --allow-provisional     # 6 生成开发划分
```

**步骤顺序不可调换。** 第 2 步会重写 `results/gv02_03_v2/family/acceptance.json` 和
`unresolved_negative_controls.tsv`；第 3 步读后者，第 5 步把前者按 SHA-256 钉进
`release_manifest.json`，第 6 步再把 `release_manifest.json` 钉进 `split_manifest.json`。
若把第 5、6 步提前到第 2、4 步之前，清单哈希立即过期，`tests/test_member_b_handoff.py`
会失败且 `--allow-provisional` 会拒绝运行。第 4 步必须重跑，否则 `acceptance.json`
不含裁决结论，`reviewed_conflict_count` 和 `scientific_acceptance` 会退回 `requires_review`。

所有 JSON 由 `write_json` 以 LF 写出（与 `write_fasta`、`write_tsv` 一致），
因此同一输入的哈希在 Windows 与 WSL 下逐字节相同，不会因检出平台漂移。

复现检查：在一次完整运行结束后执行
`python -B experiments/check_family_reproducibility.py snapshot`，再完整运行，最后执行
`python -B experiments/check_family_reproducibility.py compare`。
比较最终科学表、FASTA、HMM 和 summary，不比较含时间的工具日志。

本机可使用已安装工具的 WSL：

```powershell
wsl -d Ubuntu-20.04 --cd /mnt/d/endless/机器学习综合实践/ML-GV02-03 -- python3 -B experiments/run_family_classification.py
```

新候选使用 `--candidates data/generated/实际批次/candidates.fasta`。候选路径必须存在、非空且 ID 唯一。该入口重新生成同一输出目录的分类结果，因此不同批次应复制配置并设置独立 reference/profile/result 输出目录；默认候选仅是现有 T05 的 200 条，不代表尚未生成的新 VAE 候选。

交付文件：

| 内容 | 路径 |
| --- | --- |
| 高可信 GvpA、GvpJ/其他 Gvp 负对照及待定集 | `data/processed/gv02_03_v2/reference/*.fasta` |
| 原始记录逐条审计 | `data/processed/gv02_03_v2/reference/reference_classification.tsv` |
| 同源簇、折号、模型训练角色 | `data/processed/gv02_03_v2/reference/seed_manifest.tsv` |
| 审校来源、版本、证据代码 | `data/processed/gv02_03_v2/reference/reviewed_source_manifest.tsv` |
| 候选分类与拒判理由 | `data/processed/gv02_03_v2/candidate_family_classification.tsv` |
| GvpA/GvpJ 及其他亚家族 HMM、MSA、种子 | `models/family_profiles/` |
| 模型来源清单与每折训练/核验 ID | `models/family_profiles/profile_manifest.json` |
| 原始 HMMER/BLAST 证据、对照、汇总 | `results/gv02_03_v2/family/` |
| 最终验收检查结果 | `results/gv02_03_v2/family/acceptance.json` |

## 来源和 PF00741 的边界

原始五个本地 FASTA 共 8099 条记录，文件名和标题仅作为注释来源，不直接证明 GvpA 身份。重复序列按 SHA-256 合并，但每条原始记录保留在审计表。空序列、非标准残基、片段、预测/救回/Probable 注释及冲突注释不直接作种子。

补充来源为 UniProtKB reviewed（Swiss-Prot）GvpA/GvpJ 查询的 39 条记录。下载快照为 `references/uniprot_reviewed_gvpa_gvpj.json`，查询条件为 `reviewed:true AND (gene_exact:gvpA OR gene_exact:gvpJ)`，接口为 `https://rest.uniprot.org/uniprotkb/stream`，格式 JSON，核查日期 2026-09-22。推荐蛋白名称用于明确区分 A/J，保留片段标记和 ECO 证据代码。**人工审校不等于每条都经过功能实验。**

InterPro 官方 API `https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00741/` 的本地快照 `references/PF00741.interpro.json` 明确写道：

> This family includes gas vesicle proteins such as GvpJ.

该记录整合至 IPR000638，快照的 `subfamilies` 计数为 0。这只能证明该入口没有提供可直接使用的 A/J 子家族，不证明所有其他数据库都不存在专一模型。本项目使用已有官方 PF00741.24（Gas_vesicle，39 位点，序列及域 GA 均为 25 bits）作广义家族证据，自建 A/J HMM 作进一步竞争。未把本地旧模型称为最新 Pfam 版本。

## 参考核验方法

1. 审计标题中的明确家族标记、完整性、长度、字符与序列/登录号冲突。GvpA 种子长 50–180 aa，其他家族 40–1000 aa，属于本研究范围，不是通用生物学边界。
2. A/J 种子必须通过 PF00741 序列分数和域分数各 ≥25 bits，模型覆盖 ≥0.8。其他 Gvp 不要求 PF00741，因为 GvpN 等不属于这个域。
3. 所有合格序列使用 CD-HIT 70% identity、双向覆盖 80% 聚类。整个簇分入同一折，共三折。CD-HIT 是启发式聚类，簇隔离不保证任意跨簇序列都低于 70% 相似性。
4. A/J 的 HMM 只由完整、非 Probable 的人工审校记录作为锚点训练；其他家族用合格本地注释种子。一家族在每个簇取一条代表，降低近重复权重。每折通过 MAFFT MSA 和 hmmbuild 建模，留出簇绝不参与该折模型及竞争比对参考库。
5. 每条天然种子获得折外 HMM 分类和折外竞争性序列比对。只有原始明确家族注释与折外分类一致，才进入相应计算支持参考集。不一致记录保留审计理由，不自动改写其数据库注释。
6. 最终 A/J HMM 仍使用审校锚点，避免用模型预测自我扩充训练集。`high_confidence_gvpa.fasta` 表示当前规则下的计算高可信参考，不是实验确认集。已有生成数据 train/validation/test 不会自动替换，成员 B 应基于最终参考重新做同源簇隔离划分。

## 候选分类规则

HMMER 使用独立家族模型、`--max`、序列/域报告阈值 0 bits。只有模型覆盖及查询覆盖各 ≥0.6 的域能获选。获选家族域分数 ≥25 bits，第一、第二家族分差 ≥10 bits；局部短竞争命中仍可阻止过度确定的判定。没有可报告的第二家族时，分差留空，不伪造零分或无限分差。缺少 A/J 模型属于鉴定条件不完整，必须拒判。

HMM 高分仍可能混淆近缘 A/J。因此另用 BLASTP 竞争参考证据（SEG 开启、E-value ≤0.001、单 HSP 双向覆盖 ≥0.70；不拼接多个 HSP 重复计算覆盖）。有效 identity = 相同残基数 / max(查询长度,参考长度)，要求 ≥0.40，与次佳家族差值 ≥0.05。HMM 与序列证据必须同意；不同意或证据不足则标记歧义/不支持。比对参考库仍包含合格但未被 HMM 升级的负对照，避免把难负例过滤掉后虚报鉴定能力。

输出状态为 `supported_gvpa`、`gvpa_gvpj_ambiguous`、`other_gvp`、`unsupported`。GvpJ 属于 `other_gvp`，具体家族见 `best_family`。全部候选均有输出，包括无命中和非法序列。`supported_gvpa` 仅是家族门槛，不等于 QC 通过或具备气泡功能。

数值阈值是探索性工程规则，未以当前候选调参；本地控制数据参与了方法开发，因此报告应称内部核验，不应声称独立测试准确率。

## 正负对照和验收原则

保存全部天然折外对照，以及打乱、截短和完全重复对照。报告同时给出人工审校锚点和本地名义注释的混淆情况。完全重复对照验证分类一致性，不是独立样本，也不用于提高样本量或计算泛化准确率。

`verify_family_deliverables.py` 检查哈希、ID 与顺序、参考互斥、FASTA/TSV 一致、同源簇泄漏、重复一致性和正负对照。验收严格区分：

- 工程检查通过：文件和计算证据一致、必要对照符合检查规则。
- 科学验收通过：还要求名义非 GvpA 对照未被误纳入 GvpA；否则标记 `requires_review` 并输出具体 ID，不以隐藏难例或调高阈值掩盖失败。
- 外部准确率、全 Pfam 非目标域排查、结构和湿实验功能不在本次证据范围内。

## 结果

本轮完整运行审计 8138 条原始记录，合并为 5108 条独特序列；3063 条进入种子核验，形成 855 个 CD-HIT 簇。计算支持 GvpA 为 681 条，GvpJ 负对照 144 条，其他 Gvp 负对照 1039 条。200 条现有候选中：GvpA 支持 11 条，A/J 歧义 7 条，其他 Gvp 99 条，不支持 83 条。

**科学复核已完成，结果带有显式限定。** 1860 条名义非 GvpA 天然对照中原先 14 条被判为 GvpA；固定 NCBI 证据逐条复核后，6 条确认是本地过时的 GvpJ 注释（当前 NCBI 为 GvpA），8 条保留为 GvpA/GvpJ 操作性歧义，不再作为已确认的 GvpA。复核记录见 `results/gv02_03_v2/family/adjudication/conflict_adjudication.tsv`。

人工审校的合格去重留出记录中，14 条 GvpA 获支持、7 条证据不足；6 条 GvpJ 被正确识别为其他 Gvp、1 条证据不足，没有 GvpJ 被判为 GvpA。每类 113 条的打乱和截短对照均无 GvpA 支持，113 条完全重复对照与原序列分类一致。审校记录数量有限，且只是开发阶段内部核验，不足以证明全面泛化能力。

以 `summary.json`、`acceptance.json` 和 `control_classification.tsv` 为机器可读依据。两次完整运行输出哈希一致时，才可记录确定性复现通过。旧基线最初的 62 条误判和中间仅改用审校 HMM 的 130 条误判说明：单一 HMM 分差无法解决本数据的近缘分类问题。竞争序列证据将冲突集中到 14 条，现已逐条完成外部证据复核，并保留 8 条操作性歧义。

自动化回归：本机关闭无关 pytest 自动加载插件后运行 `python -m pytest tests -q -p no:cacheprovider`，79 项全部通过（其中分类/参考 22 项）。旧产物哈希最初因 Windows CRLF 检出失败；仅对字节归一化后与 Git 原始文件完全相同的已跟踪文本恢复 LF，并新增 `.gitattributes` 固定换行，未改写历史期望哈希。

最终工程验收：`acceptance.json` 的 12 项工程检查全部通过；两轮完整计算的 22 个最终科学产物逐文件 SHA-256 相同，见 `reproducibility_check.json`。两轮间只恢复了输入文本的原始 LF 换行，来源清单字节哈希因此不同，不混入科学结果复现比较。当前科学验收状态为 `pass_internal_controls_with_adjudicated_annotation_conflicts`；这表示冲突均已逐条裁决并保留证据，不表示独立外部准确率或实验功能已经确认。
