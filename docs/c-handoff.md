# C：B 真实 VAE 候选交接流程

本适配器读取 B 的冻结交接目录，不合并分支、不改写 B 文件、不重新分配家族标签。
输入版本为 `origin/feat/sequence-vae` 的 `10d9da8`（`member_a_v1_seed42`）：
1,000 条候选，346 条参考；主池 147、模糊池 47、排除池 806。
参考集仍属 `member_a_v1_conservative_provisional`，不是最终科学定稿。

## 运行

安装项目依赖后，在仓库根目录执行：

```powershell
.venv/Scripts/python.exe experiments/run_c_handoff.py --handoff-root .venv/b-handoff-10d9da8 --output-dir results/gv02_03_c_handoff_member_a_v1
```

`--handoff-root` 必须指向包含 B 交接文件的本地 checkout 或解压目录；这里的 `.venv` 路径
只是本机只读快照示例，不随 Git 发布。跨平台可用 `python` 替代解释器路径。
建议从提交 `10d9da8` 用 `git archive` 导出快照，保留冻结文件原始字节；Windows 自动 CRLF 转换会
改变 SHA256，此时应重新取得原始字节，不能为通过审计而修改冻结哈希。
`--config` 默认是 `configs/c_handoff_member_a_v1.json`。其他批次需独立配置及明确的新冻结输入。
输出目录必须是新的或空的，且在交接目录之外。运行时不需要 B 的 checkpoint、GPU 或重训。

## 检查与输出

- 对配置中 16 个来源文件核对 SHA256；按 ID 连接 FASTA/元数据/QC/家族/候选池表。
- 核对序列哈希、长度、EOS/上限互斥、停止原因、raw_token_length，以及批次数量。
  EOS 的 raw_token_length 为序列长度加 1，cap 为序列长度；只有 cap 标志才产生上限警告。
  这是导出记录的一致性检查，未重放模型或独立检查原始 token 流。
- 采用 B 的同一组基础 QC 参数重算并对比通过状态、原因；再增加警告模式。
  连续疏水残基至少 18、同聚物至少 6、完整串联重复至少 12 且至少 3 拷贝（原始 motif 长度 2–6）。
  警告不会改成硬性淘汰；`qc_failures` 与 `qc_warnings` 分列。
- 验证三个候选池构成无遗漏、无重复的分区，池内序列与原始候选一致；保持 B 分类及 FASTA 顺序。
- 全部 1,000 条候选分别与 346 条参考作全局 BLOSUM62 比对，保留全部最近距离并列 ID。
  结果文件 `trusted_similarity.tsv` 的名称沿用交接约定，不表示参考集已完成科学验证。
- 主池和模糊池分别提供 `*_distance.npy`、`*_distance_ids.json`、`*_pairwise.tsv`；不混合建矩阵。
  未解析距离保留 NaN（TSV 空白），不能补成 1，也不能当成多样性奖励。
- `candidate_audit.tsv` 合并 C QC、原家族/池标签及带 `global_` 前缀的比对证据。
  `manual_review.tsv` 收录 QC 失败、警告、非主池、最近参考未解析或并列任一情况。
- `qc_similarity_summary.json`、`report.md` 汇总结果；`qc_similarity_manifest.json` 记录输入/输出/
  实现文件哈希、有效参数及运行库版本。运行前后哈希不一致则拒绝发布。

可靠比对要求：双向覆盖 ≥0.8、identity ≥0.2、配对残基 ≥20、得分 >0；
distance = 1 − 相同残基数/max(两序列长度)。全局结果与 B 的局部比对证据不可混用。
单次最近参考比较上限 500,000 对；进度分块不会绕过总量限制。

## 尚未完成

额外结构域扫描及跨膜拓扑检查 **未开展**。连续疏水片段是启发式警告，不等同于跨膜预测。
QC 通过不代表家族/功能已验证，`qc_pass_extended.fasta` 不是主池。
本次不输出最终候选推荐，也不声明完整 C 生物学验证已完成。
