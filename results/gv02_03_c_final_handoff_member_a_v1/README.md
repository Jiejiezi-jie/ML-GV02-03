# C 统一复核索引（member A v1）

本目录将同一冻结 B→C 批次的主 QC/相似度、Pfam 域架构和跨膜预测逐条合并，**只增加待复核证据，不重新分池或判定功能**。完整字段解释、统计及下游使用边界见 [C→下游交接说明](../../docs/c-final-handoff.md)。

- `candidate_c_review.tsv`：1,000 条候选的原始 C 审计字段，加上域/跨膜状态和 `c_review_reasons`、`c_manual_review`。
- `manual_review_c.tsv`：`c_manual_review=True` 的 915 条子集，不是 915+37 条。
- `summary.json`：按池和告警来源计数；`manifest.json`：三个来源清单及本目录机器产物的 SHA-256。

30 条多段 PF00741 警告、7 条跨膜警告都已在原复核子集；395 条无显著 Pfam 命中维持“未解析”，不是阴性证明。原始主报告仍保留其运行当时“未检查”两项的历史记录，本目录是后续增量汇总。
