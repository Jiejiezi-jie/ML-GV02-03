# C：最近可靠天然参考匹配

该模块复用 `similarity.py` 的全局比对和可靠性门槛，完整比较每条候选与指定参考库，
输出每条候选的最近可靠参考。不改变 QC 资格、原有评分或候选间距离矩阵。

## 运行

```bash
python experiments/run_nearest_reference.py --config configs/nearest_reference.yaml --output-dir results/my_nearest_reference
```

路径相对于仓库根目录，也接受绝对路径。输出目录必须不存在或为空。
发布结果时会独占目录写入权并再次检查目录内容，避免本命令的并发运行互相覆盖。
若进程被强制终止，可能留下 `.nearest_reference.lock` 或不完整结果；请选择新的输出目录，
不要将未完成的输出当作有效运行。生成期间也不要用其他程序修改输入、代码或输出目录。
需要项目现有的 Biopython、NumPy 和 PyYAML，不依赖 BLAST/HMMER 等外部工具。
默认配置使用旧 T05 的 200 条混合候选和 1,224 条暂定 nominal-GvpA 参考。
两者均未由本模块验证家族身份，因此只能作为开发验证，不能直接认定候选属于 GvpA。
默认包含全部候选，不隐式筛除 QC 失败序列；需要通过 QC 的分析池时，应显式换输入 FASTA。

## 最近的定义

沿用全局 BLOSUM62 比对，内部和末端 gap 均采用 open=-10、extend=-0.5。
比对中同时具有残基的列数为 P，完全相同的残基对数为 M。

- `identity = M/P`；query/target coverage 分别为 P/Lq 和 P/Lr。
- `effective_identity = M/max(Lq,Lr)`，`distance = 1 - effective_identity`。
- 默认可靠性要求：双向覆盖率 ≥ 0.80、identity ≥ 0.20、P ≥ 20，且比对得分 **> 0**。
- 仅在所有门槛均通过的参考中选最小距离，不按最高原始 identity 或最高比对得分排序。
- 全部不通过时输出 `unresolved`，不强行选择“最像”的不可靠参考，不填最大距离 1。
- 对距离并列，用整数比值 M/max(Lq,Lr) 精确比较，不使用浮点容差；保留所有并列参考 ID，
  按 ID 字典序排序，首个作为代表。**表内的单组指标属于这个代表，不是所有并列参考的平均。**
- 相同序列但不同 ID 的参考全部保留；候选和参考各自要求 ID 唯一，但允许跨库同名。

比对本身的并列路径沿用现有策略：序列按字典序固定方向，取第一个最优比对。
“可靠”仅指通过这组工程门槛，不是同源关系、结构或功能的证明。

## 输出

`nearest_reference.tsv` 保持候选输入顺序，每条一行：

- `sequence_id`、`query_length`。
- `closest_reference_id`：字典序首个最近可靠参考，无可靠结果则为空。
- `closest_reference_ids`：全部最近并列 ID 的 JSON 数组；无可靠结果为 `[]`。
- `closest_reference_tie_count`：最近参考个数，单一匹配为 1，无匹配为 0。
- `identity`、`query_coverage`、`target_coverage`、`distance` 及代表比对的原始计数与得分。
- `distance_status`、`distance_reason`：无可靠匹配时为 `unresolved/no_reliable_reference`，
  此时代表参考的长度、指标、得分和距离均为空，而不是零。
- `reference_count` 和 `resolved_reference_count`：比较参考数与通过门槛的参考数。

另有 `summary.json` 和 `manifest.json`：记录候选/参考/比较数、匹配/未解析/并列数，
输入与输出 SHA-256、实际参数、版本和实现代码哈希。
本版不保存全部候选—参考比对表；可由已记录输入、参数及实现重新计算。
同一输入路径、运行环境和配置下应逐字节复现；manifest 含绝对路径与环境版本，
不承诺跨机器字节一致。代码哈希基于实际文件字节，换行转换会改变哈希。
输入和代码哈希在计算前后核验，若发生变化则拒绝发布结果。

## 资源与边界

每个库须非空，蛋白序列须为标准 20 种氨基酸，大小写统一为大写。
默认单序列最长 2,000、候选数 × 参考数最多 500,000；超预算在比对前报错。
逐候选扫描参考，不保留全部比对记录。扩大数据量前应评估运行时间，
未来可引入快速检索，但目前完整搜索不依赖预筛选命中。

不自动加入旧评分流水线，不将无可靠匹配解释为高新颖性，不执行家族归属判定。
正式分析需替换为 A 确认的参考库和相应候选池。

## 本地开发验证结果

`results/gv02_03_nearest_reference_dev/` 保存了 200 × 1,224 = 244,800 次比较后的结果：
165 条候选存在可靠参考、35 条未解析、55 条存在并列最近参考、0 条完全一致。
这 165 条是匹配门槛通过数，不是 QC 通过数，也不是确认具有功能的候选数。

两次完整运行的 `nearest_reference.tsv` 和 `summary.json` 字节一致；两次运行之间
只补了并发输出保护和源文件变化检测，因此 manifest 的实现哈希不同，不声称整个
运行目录字节一致。最终运行的输入及代码哈希均核对通过。
测试还对结果中列出的全部最近参考重新比对，检查距离及代表参考的指标一致性。
