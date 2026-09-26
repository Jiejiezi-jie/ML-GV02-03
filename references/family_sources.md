# 家族鉴定外部来源

核查日期：2026-09-22。外部 JSON 已固定在本目录，后续运行不联网刷新。

| 来源 | 请求与快照 | 用途 |
| --- | --- | --- |
| InterPro / Pfam | https://www.ebi.ac.uk/interpro/api/entry/pfam/PF00741/ → `PF00741.interpro.json` | 官方说明明确包含 GvpJ；不是 GvpA 专一域 |
| UniProtKB Swiss-Prot | https://rest.uniprot.org/uniprotkb/stream?query=reviewed%3Atrue+AND+%28gene_exact%3AgvpA+OR+gene_exact%3AgvpJ%29&format=json → `uniprot_reviewed_gvpa_gvpj.json` | 39 条人工审校原始记录，含片段/Probable 记录，运行时审计后选择种子 |
| 本地已缓存 PF00741.24 | `models/pfam/PF00741.hmm` | 39 位点、GA 25 bits；保留原版本，不宣称最新 |

每条 UniProt 记录的 accession、原始注释证据、文献、序列和版本均在 JSON 中。
可读摘要由流程生成至 `data/processed/gv02_03_v2/reference/reviewed_source_manifest.tsv`。
全部来源 SHA-256 由运行 manifest 记录。原有 RefSeq/本地 T05 注释未逐条在线重取，
因此模型与注释冲突仍需人工复核，不能据此断言数据库注释有误。
