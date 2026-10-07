#!/usr/bin/env python3
"""Build the final report and figures from hash-verified evaluation artifacts."""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gv_eval.d_evaluation import METRICS
from gv_eval.io import sha256_file, write_json
from gv_eval.reproducibility import verify_hashes
from gv_eval.selection import pareto_ranking

NAMES = ["约束满足度", "保守性", "新颖性", "候选独特性"]
LABELS = ["Constraint", "Conservation", "Novelty", "Uniqueness"]
STRATEGIES = {"weighted_sum": "等权加权", "pareto": "Pareto", "dimension_round_robin": "分维度轮转"}
COLORS = ["#2563a6", "#14918b", "#d98532"]
REPORT_ARTIFACTS = (
    "figures/correlations.png", "figures/metric_scatter_matrix.png",
    "figures/robustness.png", "figures/strategy_comparison.png",
    "pareto_front_summary.tsv", "pareto_ranking.tsv", "PROJECT_REPORT.md", "statistics.json",
)


def markdown_table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
                       "| " + " | ".join(["---"] * len(headers)) + " |",
                       *["| " + " | ".join(map(str, row)) + " |" for row in rows]])


def build(source: Path, output: Path, pareto_analysis: Path, replace_published: bool = False):
    source, output = source.resolve(), output.resolve()
    if output.exists() and any(output.iterdir()):
        if not replace_published:
            raise ValueError("Report output must be new or empty, or use --replace-published")
        previous_manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
        if set(previous_manifest["output_sha256"]) != set(REPORT_ARTIFACTS):
            raise ValueError("Existing manifest must describe exactly the generated report artifacts")
        verify_hashes(output, previous_manifest["output_sha256"])
    manifest = json.loads((source / "manifest.json").read_text(encoding="utf-8"))
    verify_hashes(source, manifest["output_sha256"])
    pareto_manifest = json.loads((pareto_analysis / "manifest.json").read_text(encoding="utf-8"))
    verify_hashes(pareto_analysis, pareto_manifest["output_sha256"])
    if pareto_manifest["source_manifest_sha256"] != sha256_file(source / "manifest.json"):
        raise ValueError("Pareto analysis and main report use different evaluation batches")
    pareto_summary = json.loads((pareto_analysis / "summary.json").read_text(encoding="utf-8"))
    pareto_comparison = pd.read_csv(pareto_analysis / "group_comparison.tsv", sep="\t")
    pareto_primary = pareto_comparison[pareto_comparison.comparison.eq("front_vs_rest")].set_index("feature")
    pareto_matched = pareto_comparison[pareto_comparison.comparison.eq("front_vs_length_matched")].set_index("feature")
    pareto_ties = pd.read_csv(pareto_analysis / "tie_match_summary.tsv", sep="\t").set_index("feature")
    read = lambda name: pd.read_csv(source / f"{name}.tsv", sep="\t")
    pool, summary = read("ranking_pool_scores"), json.loads((source / "summary.json").read_text())
    strategy, ablation = read("strategy_summary"), read("ablation_summary")
    weights, sensitivity = read("weight_robustness"), read("threshold_coverage_sensitivity")
    random = read("random_comparison")
    top20 = strategy[strategy.budget.eq(20)].set_index("strategy").loc[list(STRATEGIES)]
    pearson, spearman = (pool[METRICS].corr(method=m) for m in ("pearson", "spearman"))
    for name, actual in (("pearson", pearson), ("spearman", spearman)):
        published = pd.read_csv(source / f"correlation_{name}.tsv", sep="\t", index_col=0)
        np.testing.assert_allclose(actual, published, rtol=1e-9, atol=1e-10)
    ranked = pareto_ranking(pool, METRICS)
    frontier = ranked.loc[ranked.pareto_rank.eq(0)]
    front_ids = set(frontier.sequence_id)
    published_top = read("pareto_top20").sequence_id.tolist()
    if ranked.head(20).sequence_id.tolist() != published_top:
        raise ValueError("Stored score precision changes Pareto ordering")
    output.mkdir(parents=True, exist_ok=True)
    figures = output / "figures"
    figures.mkdir(exist_ok=True)
    ranked[["sequence_id", *METRICS, "pareto_rank", "crowding_distance"]].to_csv(
        output / "pareto_ranking.tsv", sep="\t", index=False, lineterminator="\n")
    front_stats = []
    for row in strategy.itertuples():
        selected = read(f"{row.strategy}_top{row.budget}")
        front_stats.append({"strategy": row.strategy, "budget": row.budget,
                            "pool_front_size": len(frontier),
                            "selected_from_pool_front": len(set(selected.sequence_id) & front_ids)})
    pd.DataFrame(front_stats).to_csv(output / "pareto_front_summary.tsv", sep="\t", index=False, lineterminator="\n")
    mask = pool.sequence_id.isin(front_ids).to_numpy()
    fig, axes = plt.subplots(4, 4, figsize=(10, 9))
    for i, y in enumerate(METRICS):
        for j, x in enumerate(METRICS):
            ax = axes[i, j]
            if i == j:
                ax.hist(pool[x], bins=15, color="#aac1d6", edgecolor="white")
            else:
                ax.scatter(pool.loc[~mask, x], pool.loc[~mask, y], s=9, color="#b3bec8", alpha=.6)
                ax.scatter(pool.loc[mask, x], pool.loc[mask, y], s=15, color=COLORS[1], alpha=.9)
            if i == 3:
                ax.set_xlabel(LABELS[j])
            if j == 0:
                ax.set_ylabel(LABELS[i])
            ax.tick_params(labelsize=7)
    fig.suptitle(f"Four objectives on {len(pool)} eligible candidates; teal = Pareto front ({len(frontier)})", fontsize=12)
    fig.tight_layout()
    fig.savefig(figures / "metric_scatter_matrix.png", dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8))
    for ax, corr, title in zip(axes, [pearson, spearman], ["Pearson", "Spearman"]):
        ax.imshow(corr, vmin=-1, vmax=1, cmap="RdBu_r")
        ax.set(xticks=range(4), yticks=range(4), xticklabels=LABELS, yticklabels=LABELS, title=title)
        ax.tick_params(axis="x", rotation=35)
        for i in range(4):
            for j in range(4):
                ax.text(j, i, f"{corr.iloc[i,j]:.2f}", ha="center", va="center",
                        color="white" if abs(corr.iloc[i,j]) > .6 else "#16324a")
    fig.tight_layout()
    fig.savefig(figures / "correlations.png", dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), gridspec_kw={"width_ratios": [1.25, 1]})
    x = np.arange(4)
    labels = ["Weighted", "Pareto", "Round-robin"]
    for i, (name, row) in enumerate(top20.iterrows()):
        axes[0].bar(x + (i-1)*.24, [row[f"mean_{m}"] for m in METRICS], .24, label=labels[i], color=COLORS[i])
    axes[0].set(xticks=x, xticklabels=LABELS, ylim=(0, 1.08), ylabel="Mean objective score", title="Top-20 quality")
    axes[0].legend(fontsize=8)
    for i, (_, row) in enumerate(top20.iterrows()):
        lo, hi, mean = (row[f"subset_diversity_{m}"] for m in ["lower_bound", "upper_bound", "observed_mean"])
        axes[1].hlines(i, lo, hi, color=COLORS[i], linewidth=5, alpha=.45)
        axes[1].plot(mean, i, "o", color=COLORS[i])
        axes[1].text(hi + .015, i, f"{row.subset_diversity_resolved_fraction:.1%} pairs", va="center", fontsize=9)
    axes[1].set(yticks=range(3), yticklabels=labels, xlim=(0, .85), ylim=(-.6, 2.6),
                xlabel="Diversity (distance)", title="Observed mean and unknown-pair bounds")
    axes[1].invert_yaxis()
    fig.tight_layout()
    fig.savefig(figures / "strategy_comparison.png", dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5))
    axes[0].hist(weights.top_k_jaccard, bins=np.linspace(.5, 1.01, 12), color=COLORS[0], edgecolor="white")
    axes[0].set(xlabel="Top-20 Jaccard with equal weights", ylabel="Count", title="200 weight perturbations (+/-20%)")
    threshold = sensitivity[sensitivity.score_factor.eq(1)].pivot(index="coverage_offset", columns="minimum_resolved_fraction", values="jaccard")
    counts = sensitivity[sensitivity.score_factor.eq(1)].pivot(index="coverage_offset", columns="minimum_resolved_fraction", values="ranking_count")
    axes[1].imshow(threshold, vmin=0, vmax=1, cmap="Blues", aspect="auto")
    axes[1].set(xticks=range(3), xticklabels=[f"{x:.0%}" for x in threshold.columns], yticks=range(3),
                yticklabels=[f"{min(.95+x,1):.2f}" for x in threshold.index],
                xlabel="Required reliable-pair fraction", ylabel="Domain coverage threshold", title="Threshold sensitivity (score cutoff 25)")
    for i in range(3):
        for j in range(3):
            axes[1].text(j, i, f"J={threshold.iloc[i,j]:.2f}\nn={counts.iloc[i,j]}", ha="center", va="center", color="white")
    fig.tight_layout()
    fig.savefig(figures / "robustness.png", dpi=170)
    plt.close(fig)
    stats = {
        "source_manifest_sha256": sha256_file(source / "manifest.json"),
        "candidate_count": summary["funnel"]["input_candidate_count"],
        "eligible_count": summary["funnel"]["eligible_candidate_count"],
        "ranking_count": len(pool), "pareto_front_size": len(frontier),
        "pareto_layer_count": int(ranked.pareto_rank.max()) + 1,
        "weight_jaccard": {"mean": float(weights.top_k_jaccard.mean()), "min": float(weights.top_k_jaccard.min()), "max": float(weights.top_k_jaccard.max())},
        "weight_rank_spearman_mean": float(weights.rank_spearman.mean()),
        "threshold_jaccard": {"min": float(sensitivity.jaccard.min()), "max": float(sensitivity.jaccard.max())},
        "threshold_ranking_count": {"min": int(sensitivity.ranking_count.min()), "max": int(sensitivity.ranking_count.max())},
        "top20": top20.reset_index().to_dict(orient="records"), "front_by_strategy": front_stats,
        "ablation": ablation.fillna("").to_dict(orient="records"),
        "pearson": pearson.to_dict(), "spearman": spearman.to_dict(),
        "pareto_analysis_manifest_sha256": sha256_file(pareto_analysis / "manifest.json"),
        "optional_experiments": {"5": "not_done", "6": "completed", "7": "completed"},
    }
    write_json(output / "statistics.json", stats)
    quality_rows = [[STRATEGIES[name], *[f"{row['mean_'+m]:.4f}" for m in METRICS],
                     f"{row.subset_diversity_observed_mean:.4f}", f"{row.subset_diversity_resolved_fraction:.1%}",
                     f"[{row.subset_diversity_lower_bound:.4f}, {row.subset_diversity_upper_bound:.4f}]"] for name, row in top20.iterrows()]
    table_quality = markdown_table(["策略", *NAMES, "已解析对多样性", "配对覆盖率", "全部配对均值界限"], quality_rows)
    table_ablation = markdown_table(["移除维度", "名单 Jaccard", "约束均值", "新颖性均值"],
        [[NAMES[METRICS.index(row.removed_metric)], f"{row.jaccard_with_full:.4f}", f"{row.mean_constraint_score:.4f}", f"{row.mean_novelty_score:.4f}"] for row in ablation.itertuples() if isinstance(row.removed_metric, str)])
    table_front = markdown_table(["策略", "Top-10 前沿成员", "Top-20 前沿成员", "Top-50 前沿成员"],
        [[STRATEGIES[name], *[next(r['selected_from_pool_front'] for r in front_stats if r['strategy']==name and r['budget']==k) for k in (10,20,50)]] for name in STRATEGIES])
    pair_rows = []
    for i in range(4):
        for j in range(i+1, 4):
            pair_rows.append([NAMES[i] + " / " + NAMES[j], f"{pearson.iloc[i,j]:.4f}", f"{spearman.iloc[i,j]:.4f}"])
    random_rows = []
    for name in STRATEGIES:
        for metric in ["mean_constraint_score", "mean_novelty_score", "subset_diversity_lower_bound", "subset_diversity_resolved_fraction"]:
            r = random[(random.strategy==name)&(random.metric==metric)].iloc[0]
            label = {"mean_constraint_score":"约束", "mean_novelty_score":"新颖性", "subset_diversity_lower_bound":"多样性下界", "subset_diversity_resolved_fraction":"配对覆盖率"}[metric]
            random_rows.append([STRATEGIES[name], label, f"{r.observed:.4f}", f"{r.random_mean:.4f}", f"{r.empirical_p_random_ge_observed:.4f}"])
    report = f'''# GvpA 候选序列的多目标评价与筛选

## 摘要

本项目为 GvpA 蛋白候选构建可复现的评价与筛选系统，研究四个目标的冲突、聚合策略差异和参数敏感性。发布批次包含 1,000 条序列 VAE 候选。质量、家族、结构域联合门槛保留 136 条，其中 114 条具备足够的四维证据，进入相同候选池上的三策略比较。系统导出 Top-10、Top-20、Top-50 名单，完成逐维消融、200 次权重扰动、27 组门槛组合和 1,000 次随机筛选对照。

在该批次中，保守性与候选独特性冲突最明显，新颖性与候选独特性高度相关。等权加权保留较高的约束分，Pareto 更偏向探索新颖候选。可靠配对比例的门槛比域分数的局部扰动更影响名单。结果支持按验证目标选择策略，并同时检查证据覆盖与人工复核标记。所有结论均为固定暂定参考集上的计算评价。

## 1. 研究问题与任务范围

- RQ1：约束满足度、保守性、新颖性和多样性的冲突结构如何？
- RQ2：在同一候选池、同一预算下，三种聚合策略如何改变名单和质量？
- RQ3：权重与门槛的局部扰动是否导致明显变化？

课程原文见 [项目要求第 20–22 页](../docs/course/project-pool.pdf)。GV02-03 要求构建评价系统，不要求训练新模型。仓库中的 GRU-VAE 是候选来源扩展；评价和筛选可直接使用已冻结序列。选做采用实验 6（前沿目标取舍与序列特征）和实验 7（同池随机筛选对照），共两项；实验 5 下游性质预测未做。

## 2. 数据与筛选流程

346 条保守参考按同源簇拆分为 277/33/36 条训练、验证和测试序列，31 个簇互不跨集合。生成批次固定种子 42，包含 1,000 条唯一序列，其中 972 条自然结束，28 条达到生成长度上限。原训练权重尚未入库，训练日志和生成元数据随数据发布。

| 阶段 | 数量 | 解释 |
| --- | ---: | --- |
| 原始候选 / 基础 QC 通过 | 1,000 / 1,000 | 格式、长度及基础序列检查 |
| GvpA 家族支持 | 147 | 另外 47 条歧义，806 条排除 |
| QC、家族、域联合合格 | 136 | PF00741 域分数 ≥25 bits、模型覆盖 ≥0.95 |
| 四维证据充分 | 114 | 22 条评分证据不足，仍保留审计记录 |

PF00741 也覆盖部分 GvpJ，家族身份采用竞争性证据单独判定。参考发布状态仍为 `provisional_not_scientifically_final`。原始来源、人工争议复核、家族模型和逐条输出保留在数据及证据目录中。

## 3. 方法

### 3.1 四个目标

1. **约束满足度**：PF00741 域分数相对于参考分布的经验百分位，结合模型覆盖因子。域门槛先决定资格，连续分数再用于排序。
2. **保守性**：固定参考 MSA 的 23 个统计保守位点上，候选与共识残基的一致比例。位点来自参考序列，候选不参与位点定义。
3. **新颖性**：候选到可靠天然参考的最小全局距离。全局比对使用 BLOSUM62；有效一致率为一致残基数除以较长序列长度，距离为 1 减该一致率。低覆盖、低一致率或低分的比对标记为未解析。
4. **候选独特性**：候选到 136 条联合门槛合格序列中其他条目的已解析距离均值。进入排名需至少 50% 的配对可靠；所有 114 条排名候选使用同一个参照池。

候选独特性用于逐条排序，集合多样性另外计算入选集合内部的两两距离。原 147 条主池共有 10,731 对，其中 3,377 对未解析。缺失距离保留 NaN，不能用最大距离奖励未解析候选。报告同时给出已解析对的条件均值、有效配对比例以及所有未知距离位于 [0,1] 时的均值上下界；该界限不是置信区间。

### 3.2 策略与实验设置

四维分数均位于 [0,1]，越大越好。等权加权使用四个 0.25 权重；Pareto 按非支配层、拥挤距离和稳定 ID 排序；分维度轮转按各维度排序依次纳入尚未入选候选。预算为 K=10/20/50，主分析 K=20。家族、QC、域门槛对三种策略完全一致，名单不足时明确报错或标记，不从其他池补足。

逐维消融仅移除排序目标，保留硬门槛。权重在各基准权重附近独立扰动 ±20% 后重新归一化，重复 200 次。域分数采用 0.8/1.0/1.2 倍阈值，域覆盖偏移 -0.1/0/+0.1（截断到 1），可靠配对要求为 25%/50%/75%，共 27 组。随机对照在同一 114 条候选池无放回抽取 20 条，共 1,000 次。随机种子均为 42。

25 bits 来自 Pfam 模型 GA，域覆盖由参考分布与预设上下限共同校准；50% 配对比例是显式工程门槛，其合理性通过敏感性分析呈现。

## 4. RQ1：目标冲突

{markdown_table(['维度对', 'Pearson', 'Spearman'], pair_rows)}

![相关性](figures/correlations.png)

![四维散点矩阵](figures/metric_scatter_matrix.png)

当前池中，保守性与候选独特性是负相关最强的一对，Pearson 为 {pearson.iloc[1,3]:.4f}。新颖性与候选独特性的 Pearson 为 {pearson.iloc[2,3]:.4f}，说明两个目标在当前批次中有较强重叠，但定义与参照对象仍不同。相关性只描述经过硬门槛筛选后的 114 条，不推广到全部生成序列或天然蛋白。

## 5. RQ2：策略比较与 Pareto 前沿

{table_quality}

![策略质量和集合多样性](figures/strategy_comparison.png)

等权加权在 Top-20 的约束和保守性均值最高；Pareto 的新颖性和候选独特性均值最高。Pareto 的条件多样性均值也更高，但配对覆盖率更低，且三种策略的未知距离界限重叠，因此不能据此断言完整集合多样性严格更高。

114 条候选的第一非支配前沿包含 **{len(frontier)} 条**，共有 {stats['pareto_layer_count']} 个非支配层。前沿成员指在这四个目标上不被同池其他候选支配的序列，不能解释为功能最优蛋白。各策略从整个候选池第一前沿选入的数量如下。

{table_front}

完整九份名单及重叠度见 [筛选结果目录](../results/gv02_03_v2/d_evaluation_member_a_v1/) 和其中的 `strategy_overlap.tsv`；[完整 Pareto 层级](pareto_ranking.tsv) 保留每条候选的四维分数。

## 6. 消融与 RQ3：鲁棒性

{table_ablation}

移除约束分使 Top-20 名单变化最大。移除保守性影响最小，与本批次保守性已接近饱和有关；这不能证明该维度在其他数据或取消硬门槛后冗余。新颖性和独特性虽高度相关，删除任一项仍会改变名单。

200 次权重扰动的 Jaccard 均值为 **{stats['weight_jaccard']['mean']:.4f}**，范围为 {stats['weight_jaccard']['min']:.4f}–{stats['weight_jaccard']['max']:.4f}；全池排序 Spearman 均值为 {stats['weight_rank_spearman_mean']:.4f}。27 组门槛实验的 Jaccard 为 {sensitivity.jaccard.min():.4f}–{sensitivity.jaccard.max():.4f}，可排名候选为 {sensitivity.ranking_count.min()}–{sensitivity.ranking_count.max()} 条。

![参数敏感性](figures/robustness.png)

固定 50% 可靠配对比例时，本次域分数和域覆盖的全部扰动都保持加权 Top-20 名单不变；改为 25% 或 75% 时名单明显变化。门槛变化会重建候选独特性的参照池，因而同时影响资格和分数。提高稳定性的办法是固定参考与比对规则、预先声明覆盖门槛、同时报告多个门槛下的名单及入选频率，并优先复核未知配对较多的候选。

## 7. 选做：与随机筛选对比

{markdown_table(['策略', '指标', '策略值', '随机均值', '随机不低于策略的经验比例'], random_rows)}

三种策略均提高本批次的新颖性和候选独特性，但 Pareto 与轮转的约束均值低于随机均值，配对覆盖率也低于随机。因此多目标筛选并非每一维都优于随机。经验比例采用加一修正，仅描述固定候选池上的随机抽样对照，未做多重比较校正，不代表生物学功能显著改善。

## 8. 选做实验 6：前沿取舍与序列共同特征

在固定的 114 条候选上重算四维前沿，对全部 10 条成员列出原始分数、同池百分位和复核证据，并与 104 条非前沿比较。分析前固定长度、组成熵、最高单一残基比例、DEKR 比例、G/P 比例及最长同聚物六项描述，完整报告 20 种残基的组成。另用最小总绝对长度差的一对一匹配建立 10 条对照，8 对等长，其余各差 1 个残基。

![四维前沿及逐条取舍](../results/gv02_03_v2/pareto_analysis_member_a_v1/figures/frontier_objectives.png)

`vae_candidate_000959` 的约束分为 1.0000；`vae_candidate_000061` 的新颖性和独特性分别为 0.5818 与 0.6248，但约束分为 0.2787。两者体现不同目标取舍，前沿成员不代表功能最优。

前沿长度均值为 {pareto_primary.loc['length','front_mean']:.1f} 个残基，其余为 {pareto_primary.loc['length','control_mean']:.2f}。长度匹配后，组成熵均值差为 {pareto_matched.loc['shannon_entropy_bits','mean_difference']:+.5f} bits，最高单一残基比例差为 {100*pareto_matched.loc['max_residue_fraction','mean_difference']:+.3f} 个百分点。G/P 比例的组间均值差由 {100*pareto_primary.loc['gly_pro_fraction','mean_difference']:+.3f} 变为 {100*pareto_matched.loc['gly_pro_fraction','mean_difference']:+.3f} 个百分点，说明这一差异依赖对照选择。

进一步检查 {pareto_summary['length_match_sensitivity']['runs']} 组等价最优匹配，每组总长度差均为 2。组成熵差始终为正（{pareto_ties.loc['shannon_entropy_bits','min_mean_difference']:+.5f} 至 {pareto_ties.loc['shannon_entropy_bits','max_mean_difference']:+.5f} bits），最高单一残基比例差始终为负。G/P 差异有 {int(pareto_ties.loc['gly_pro_fraction','negative_count'])} 组为负、{int(pareto_ties.loc['gly_pro_fraction','positive_count'])} 组为正，因此不能把更高 G/P 比例作为稳健的前沿共同特征。这是事后增加的对照选择敏感性检查，不是独立重复或显著性检验。

前沿仍有 {pareto_summary['evidence']['front']['manual_review_counts']['c_manual_review']}/10 条需人工复核，逐候选可靠配对比例均值为 {pareto_summary['evidence']['front']['uniqueness_resolved_fraction']['mean']:.2%}，其余为 {pareto_summary['evidence']['rest']['uniqueness_resolved_fraction']['mean']:.2%}。组成特征仅作当前批次的描述，不是下游性质预测或功能验证。完整特征、长度匹配及等价匹配敏感性结果见[实验 6 报告](../results/gv02_03_v2/pareto_analysis_member_a_v1/EXPERIMENT_6_REPORT.md)。

## 9. 适用场景、复现与边界

- 验证预算有限、重视域支持时：以加权名单为起点，根据需求预先设定权重，并逐条审阅预警。
- 希望探索目标间取舍时：查看 Pareto 前沿及原始四维证据，关注新颖性提升伴随的约束下降和未知距离。
- 希望覆盖各维度高分候选时：采用分维度轮转，保持固定总预算，并检查其整体均衡性。

运行环境、命令与证据范围见 [复现说明](../docs/REPRODUCIBILITY.md)。评价入口重算评分、筛选、消融、鲁棒性和随机基线；扩展复现入口另外重算 QC 和全局比对，并核验发布的家族、Pfam 和跨膜证据。

原始 VAE checkpoint 缺失，因此原训练和原采样尚不能重放；家族模型构建和完整 Pfam/PureseqTM 扫描需要额外外部工具及数据库。本项目未验证候选的表达、结构、组装或功能，结果用于后续实验排序。课程要求中的评价工具与计算实验不依赖原始生成权重。

### 数据来源

基础结果取自 `results/gv02_03_v2/d_evaluation_member_a_v1/` 的冻结产物，其输入来自已固定的生成与全局比对提交；实验 6 取自同批次的 `pareto_analysis_member_a_v1/`。本文及图表由 `experiments/build_project_report.py` 读取经过 SHA-256 校验的文件生成。基础源清单 SHA-256：`{stats['source_manifest_sha256']}`；实验 6 清单 SHA-256：`{stats['pareto_analysis_manifest_sha256']}`。
'''
    root_link = Path(os.path.relpath(ROOT, output)).as_posix()
    for directory in ("docs", "results"):
        report = report.replace(f"(../{directory}/", f"({root_link}/{directory}/")
    (output / "PROJECT_REPORT.md").write_text(report, encoding="utf-8", newline="\n")
    write_json(output / "manifest.json", {
        "source_manifest_sha256": stats["source_manifest_sha256"],
        "pareto_analysis_manifest_sha256": stats["pareto_analysis_manifest_sha256"],
        "builder_sha256": sha256_file(Path(__file__)),
        "output_sha256": {name: sha256_file(output / name) for name in REPORT_ARTIFACTS},
    })
    print(json.dumps({"report": str(output / 'PROJECT_REPORT.md'), "pareto_front_size": len(frontier),
                      "weight_jaccard": stats['weight_jaccard']}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=ROOT / "results/gv02_03_v2/d_evaluation_member_a_v1")
    parser.add_argument("--output", type=Path, default=ROOT / "reports")
    parser.add_argument("--pareto-analysis", type=Path, default=ROOT / "results/gv02_03_v2/pareto_analysis_member_a_v1")
    parser.add_argument("--replace-published", action="store_true", help="Verify the existing report manifest, then replace only generated report artifacts")
    args = parser.parse_args()
    build(args.source, args.output, args.pareto_analysis, args.replace_published)
