# C 模块：全局相似度与距离

## 运行入口

```bash
python -B experiments/run_similarity.py --config configs/similarity.yaml
```

依赖项目已声明的 Biopython、NumPy 和 PyYAML。输入为一个独立候选池的 FASTA；
配置中的 `pool_label` 记录调用方对该池的描述，不作为家族鉴定证据。
默认配置使用原始 T05 200 条混合候选，仅用于开发验证。
正式运行应由 A/D 提供通过 QC 和家族划分的 FASTA，分别运行主分析池与歧义池。

输出目录必须为空或不存在，重复运行需配置新的目录，防止覆盖已有证据。
序列须非空、仅含标准 20 种氨基酸，ID 唯一；允许不同 ID 有相同序列。
大小写统一为大写，原始 FASTA 不修改。空池报错，单序列池输出 1×1 零矩阵和仅含表头的成对表。

## 比对和指标

采用 Biopython `PairwiseAligner` 的全局模式、BLOSUM62 和 affine gap。
默认 gap open=-10、gap extend=-0.5；末端缺口和内部缺口使用相同罚分。
选择一个最优比对，固定按序列字符串的字典序决定比对方向，避免交换输入时
由于同分比对选择不同而产生非对称距离；同分时取固定后端返回的第一个比对。
manifest 保存后端版本和选择规则。不同后端版本的同分比对可能不同。

参考实现说明：[Biopython pairwise alignment](https://biopython.org/docs/latest/Tutorial/chapter_pairwise.html)。

设两条序列长度为 Lq、Lt，比对中双方都有残基的列数为 P，完全相同的残基对数为 M：

| 字段 | 定义 |
| --- | --- |
| `alignment_score` | BLOSUM62 与 gap penalty 的总分；不是概率或 E-value |
| `alignment_length` | 完整比对列数，包括缺口列 |
| `aligned_pairs` | P，不包括任一侧是缺口的列 |
| `identical_residues` | M |
| `identity` | M/P，P=0 时为 0；取值为 0–1，不是百分数 |
| `query_coverage` | P/Lq |
| `target_coverage` | P/Lt |
| `effective_identity` | M/max(Lq,Lt) |
| `distance` | 通过可靠性规则后为 1−effective_identity，否则缺失 |

全局比对会消费全部输入字符，但覆盖率仍仅统计双方都有残基的列，
因此长末端或内部插入不会自动获得 100% 双向覆盖。

## 可靠性规则与缺失值

默认同时要求：双向覆盖率 ≥0.80、identity ≥0.20、P ≥20，以及 alignment score >0。
这些是冻结在配置中的初始工程筛查规则，尚未经 GvpA/GvpJ 正负对照校准。
`resolved` 仅表示通过这些规则，不证明同源、家族身份或功能。
阈值来自配置；非法、非有限数值和超范围参数直接报错。

失败的每个条件记录于 `distance_reason`，状态为 `unresolved`；仍保留 identity、覆盖率等原始证据。
TSV 中缺失 distance 留空，NPY 中为 NaN；不得将它补成 1，也不得直接用 `nanmean`
忽略缺失后声称得到了完整的集合多样性。D 应另行制定不完整比较的排除/报告策略。
对角线固定为 0，是自距离约定，不经过可靠性门槛。
该距离是覆盖校正的不相似度，不保证满足数学度量的三角不等式。

## 输出接口与审计

- `pairwise_similarity.tsv`：每个无序序列对一行，包含双方 ID、长度、全部比对指标和状态/原因。
- `candidate_distance.npy`：对称 float64 矩阵，禁止 pickle，缺失值为 NaN。
- `candidate_distance_ids.json`：严格保持输入 FASTA 顺序，对应矩阵行列。
- `summary.json`：候选数、成对数、resolved/unresolved 数及比例，明确未执行家族鉴定。
- `manifest.json`：输入 FASTA、配置、输出和核心源码的 SHA-256；有效参数、算法和环境版本。

所有序列对均做全局比对，因此本版没有 BLAST 无命中/Multi-HSP 聚合问题，也无需安装外部 BLAST。
这是供 C/D 使用的独立全局比对模块，未替换旧实验的 BLAST 评分逻辑。
本版尚不包含 BLAST/DIAMOND 两阶段检索、候选对天然参考的最近邻查找、半全局模式或自动家族分池。
较大集合仍需后续增加快速检索。全对全计算规模为 N(N−1)/2；默认最多 500000 对、单序列最长 2000 aa，
超过限制会在比对前报错。运行所需时间还随序列长度增长。

## 验证

2026-09-22 开发运行：旧 T05 200 条混合候选共 19900 对，6710 对通过当前门槛，
13190 对为 unresolved。两次独立运行的成对表、NPY 矩阵、ID 顺序和摘要 SHA-256 一致，
记录位于 `results/gv02_03_similarity_dev/reproducibility_check.json`。
两次配置和输出目录不同，因此不宣称两个 manifest 文件逐字节相同。
这些结果不属于正式 GvpA 同家族多样性实验。

仓库对该开发结果目录禁用 Git 文本换行转换，以保持产物字节哈希。
manifest 的源码哈希记录实际运行文件字节（包含本地换行），不同换行检出的源码哈希可能不同。

`tests/test_similarity.py` 包含手工可核算的相同序列、单替换、内部插入、长末端、
无关序列、阈值边界、反向对称性、单序列矩阵、非法输入和资源限制测试。
命令行测试两次独立运行并比较四个核心产物 SHA-256，验证输入不修改、NaN 保留及拒绝覆盖已有结果。
