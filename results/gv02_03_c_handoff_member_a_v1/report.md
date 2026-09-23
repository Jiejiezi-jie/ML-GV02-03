# C：真实 VAE 候选交接检查

开发性结果，非最终科学结论。候选 1000 条，参考 346 条。

基础 QC：1000 通过，0 失败；223 条有警告。
EOS：972；生成上限：28。仅核验导出元数据、哈希和长度一致性，未重放模型。

## 结果

候选池（沿用 B）：{"main_supported_gvpa": 147, "ambiguous_exploration": 47, "excluded": 806}。

全局最近参考：801 可解析，199 未解析，262 条存在并列。

警告计数（可重叠）：{"composition_outlier_warning": 174, "generation_cap_warning": 28, "homopolymer_warning": 46}。

人工复核清单 915 条：QC 失败、警告、非主池、未解析或并列之一。

## 使用边界

全局 BLOSUM62 比对；distance = 1 - M/max(Lq,Lr)，只对可靠比对赋值。不是 B 的局部比对结果；最近参考不能替代家族鉴定。

main / ambiguous 分别提供矩阵、顺序 ID 和逐对证据；未解析为 NaN（TSV 留空），不可补成 1。qc_pass_extended.fasta 含所有基础 QC 通过者，不等于主候选池。

**尚未检查：额外结构域、跨膜拓扑。** 连续疏水片段只作警告，不是跨膜区预测；家族标签沿用 B，没有重新验证或做最终候选推荐。
