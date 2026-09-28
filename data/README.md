# 数据

## 当前研究批次

| 目录 | 用途 |
| --- | --- |
| `raw/` | 天然序列来源与经审阅条目，保留原始 FASTA |
| `generated/sequence_vae/member_a_v1_seed42/` | 1,000 条生成序列、逐条采样元数据与生成清单 |
| `processed/gv02_03_v2/reference/` | 竞争性家族鉴定后的参考分类、负对照与种子 |
| `processed/gv02_03_v2/member_b_handoff/` | 346 条保守参考及其发布清单 |
| `processed/gv02_03_v2/member_b_handoff/development_split/` | 按同源簇隔离的 277/33/36 划分及词表 |
| `processed/gv02_03_v2/vae_member_a_v1/` | 候选 QC、家族分类、分池及原始评分证据 |

FASTA 标识与表格 `sequence_id` 一一对应，清单保存文件与序列哈希。歧义、排除和证据不足条目均保留，不能从可排名池反推所有生成候选都通过检查。

`candidates/t05/generated_gvp.fasta` 的 200 条旧候选继续作为家族鉴定历史输入保存，当前报告和筛选入口使用 1,000 条 VAE 批次。该文件属于来源依赖，不是第二套最终结果。

参考状态为 `provisional_not_scientifically_final`。PF00741 域命中不能独立证明 GvpA 身份，具体规则见[方法说明](../docs/METHODS.md)。
