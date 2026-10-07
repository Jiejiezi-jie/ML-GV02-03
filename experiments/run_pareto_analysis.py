#!/usr/bin/env python3
"""Analyze the fixed four-objective Pareto front and its sequence characteristics."""
from __future__ import annotations

import argparse
import json
import platform
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gv_eval.io import sha256_file, write_json
from gv_eval.pareto_features import AA_COLUMNS, FEATURES, METRICS, analyze_pareto_features, read_fasta_strict
from gv_eval.reproducibility import verify_hashes
from gv_eval.selection import pareto_ranking

SOURCE = ROOT / "results/gv02_03_v2/d_evaluation_member_a_v1"
OUTPUT = ROOT / "results/gv02_03_v2/pareto_analysis_member_a_v1"
NAMES = dict(zip(METRICS, ["约束满足度", "保守性", "新颖性", "候选独特性"]))
LABELS = ["Constraint", "Conservation", "Novelty", "Uniqueness"]
PALETTE = ["#168b82", "#aab6c1", "#437ab1"]


def table(headers, rows):
    return "\n".join(["| " + " | ".join(headers) + " |",
                      "| " + " | ".join(["---"] * len(headers)) + " |",
                      *["| " + " | ".join(map(str, row)) + " |" for row in rows]])


def save_figure(fig, stem, qa_scripts=None):
    """Export deterministic figures, measuring final panel geometry before export."""
    fig.canvas.draw()
    bounds = [list(ax.get_position().bounds) for ax in fig.axes if ax.get_label() != "<colorbar>"]
    write_json(stem.parent / (stem.name + ".geometry.json"), {
        "figure_inches": list(map(float, fig.get_size_inches())),
        "axes_bounds_figure_fraction": bounds,
    })
    if qa_scripts:
        sys.path.insert(0, str(qa_scripts))
        from audit_panel_alignment import require_matplotlib_panel_alignment
        qa = stem.parent.parent / "qa"
        qa.mkdir(exist_ok=True)
        require_matplotlib_panel_alignment(fig, json_out=qa / (stem.name + ".alignment.json"),
                                          tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True)
    fig.savefig(stem.with_suffix(".png"), dpi=300, facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), metadata={"Date": None})
    fig.savefig(stem.with_suffix(".pdf"), metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)


def make_figures(result, output, qa_scripts=None):
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"], "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none", "svg.hashsalt": "gv02-experiment6-v1",
                         "pdf.fonttype": 42, "legend.frameon": False})
    figures = output / "figures"
    figures.mkdir()
    all_rows = result["candidate_features"]
    front = result["front_details"].sort_values("sequence_id")
    fig, (space_ax, ax) = plt.subplots(2, 1, figsize=(9.0, 9.8))
    nonfront = all_rows.loc[~all_rows.is_front]
    for _, row in nonfront.iterrows():
        space_ax.plot(range(4), row[list(METRICS)], color="#b6c1cc", linewidth=.65, alpha=.35)
    for _, row in front.iterrows():
        space_ax.plot(range(4), row[list(METRICS)], color=PALETTE[0], linewidth=1.5, alpha=.8)
    space_ax.set(xticks=range(4), xticklabels=LABELS, ylim=(-.04, 1.04),
                 ylabel="Original objective score", title="a  Four-dimensional objective space: 10 front / 104 remaining")
    space_ax.plot([], [], color=PALETTE[0], label="First front (n=10)")
    space_ax.plot([], [], color="#aab6c1", label="Remaining pool (n=104)")
    space_ax.legend(loc="upper center", bbox_to_anchor=(.5, -.13), ncol=2, fontsize=8)
    values = front[[m + "_percentile" for m in METRICS]].to_numpy()
    im = ax.imshow(values, vmin=0, vmax=100, cmap="Blues", aspect="auto")
    ax.set(xticks=range(4), xticklabels=LABELS, yticks=range(len(front)),
           yticklabels=[s.replace("vae_candidate_", "") for s in front.sequence_id],
           ylabel="Candidate ID suffix")
    ax.set_title("b  Individual trade-offs: percentile color, raw-score labels", pad=14)
    for i in range(len(front)):
        for j, metric in enumerate(METRICS):
            ax.text(j, i, f"{front.iloc[i][metric]:.3f}", ha="center", va="center",
                    color="white" if values[i, j] >= 70 else "#17334d", fontsize=10)
    fig.subplots_adjust(left=.14, right=.84, bottom=.08, top=.95, hspace=.38)
    colorbar_ax = fig.add_axes([.87, .08, .022, ax.get_position().height], label="<colorbar>")
    colorbar = fig.colorbar(im, cax=colorbar_ax)
    colorbar.set_label("Percentile in the fixed pool (n=114)")
    save_figure(fig, figures / "frontier_objectives", qa_scripts)

    nonfront = all_rows.loc[~all_rows.is_front]
    control_ids = result["matched_pairs"].control_id
    matched = all_rows.set_index("sequence_id").loc[control_ids]
    titles = ["Sequence length", "Composition entropy", "Largest residue fraction",
              "DEKR fraction", "G/P fraction", "Longest homopolymer"]
    units = ["Residues", "Bits", "Fraction", "Fraction", "Fraction", "Residues"]
    fig, axes = plt.subplots(2, 3, figsize=(11.4, 6.8))
    # Only the horizontal display jitter is random; every measured value is retained.
    rng = np.random.default_rng(42)
    for index, (ax, feature, title, unit) in enumerate(zip(axes.flat, FEATURES, titles, units)):
        groups = [front[feature].to_numpy(), nonfront[feature].to_numpy(), matched[feature].to_numpy()]
        for i, (values, color) in enumerate(zip(groups, PALETTE)):
            ax.scatter(i + rng.uniform(-.12, .12, len(values)), values, s=17 if i != 1 else 10,
                       color=color, alpha=.85 if i != 1 else .45, edgecolors="none", zorder=2)
            q1, median, q3 = np.quantile(values, [.25, .5, .75])
            ax.plot([i, i], [q1, q3], color="#18324b", linewidth=2.5, zorder=3)
            ax.plot([i-.14, i+.14], [median, median], color="#18324b", linewidth=2, zorder=3)
        ax.set(xticks=[0, 1, 2], xticklabels=["Front\n(n=10)", "Rest\n(n=104)", "Matched\n(n=10)"],
               ylabel=unit, title=title, xlim=(-.5, 2.5))
        ax.text(-.15, 1.10, chr(97 + index), transform=ax.transAxes, fontweight="bold", fontsize=11)
        ax.margins(y=.12)
    fig.suptitle("Predefined sequence features: all candidates, medians and interquartile ranges", fontsize=12, y=.98)
    fig.subplots_adjust(left=.08, right=.98, bottom=.10, top=.88, hspace=.56, wspace=.38)
    save_figure(fig, figures / "sequence_features", qa_scripts)

    fig, ax = plt.subplots(figsize=(10, 4.6))
    x = np.arange(len(AA_COLUMNS))
    for offset, group, color, name in [(-.18, nonfront, "#718b9c", "Front - remaining 104"),
                                       (.18, matched, "#168b82", "Front - length-matched 10")]:
        difference = (front[list(AA_COLUMNS)].mean() - group[list(AA_COLUMNS)].mean()) * 100
        ax.bar(x + offset, difference, width=.34, color=color, label=name)
    ax.axhline(0, color="#51616d", linewidth=.7)
    ax.set(xticks=x, xticklabels=[c[3:] for c in AA_COLUMNS], xlabel="Amino acid (complete set of 20)",
           ylabel="Mean composition difference (percentage points)")
    ax.legend(loc="upper center", bbox_to_anchor=(.5, 1.20), ncol=2)
    fig.subplots_adjust(left=.10, right=.98, top=.80, bottom=.17)
    save_figure(fig, figures / "amino_acid_composition", qa_scripts)


def write_report(result, output, source_hash):
    front = result["front_details"].sort_values("sequence_id")
    comparison = result["group_comparison"]
    feature_names = dict(zip(FEATURES, ["长度", "组成熵（bits）", "最高单一残基比例", "DEKR 比例", "G/P 比例", "最长同聚物段"]))
    rows = []
    for row in front.itertuples():
        rows.append([row.sequence_id, *[f"{getattr(row,m):.4f}" for m in METRICS],
                     "、".join(NAMES[m] for m in row.strongest_metrics.split(";")),
                     "、".join(NAMES[m] for m in row.weakest_metrics.split(";"))])
    detail = table(["候选 ID", *NAMES.values(), "最高百分位维度", "最低百分位维度"], rows)
    # Column names follow the pure analysis module's documented API.
    primary = comparison[comparison.comparison.eq("front_vs_rest")].set_index("feature")
    matched = comparison[comparison.comparison.eq("front_vs_length_matched")].set_index("feature")
    features_table = table(["特征", "前沿中位数", "非前沿中位数", "均值差", "Cliff delta", "匹配后均值差"], [
        [feature_names[f], f"{primary.loc[f,'front_median']:.4f}", f"{primary.loc[f,'control_median']:.4f}",
         f"{primary.loc[f,'mean_difference']:.4f}", f"{primary.loc[f,'cliffs_delta']:.4f}",
         f"{matched.loc[f,'mean_difference']:.4f}"] for f in FEATURES])
    aa_table = table(["残基", "前沿均值", "非前沿均值", "主比较差（百分点）", "长度匹配差（百分点）"], [
        [a[3:], f"{primary.loc[a,'front_mean']:.4f}", f"{primary.loc[a,'control_mean']:.4f}",
         f"{100*primary.loc[a,'mean_difference']:.3f}", f"{100*matched.loc[a,'mean_difference']:.3f}"] for a in AA_COLUMNS])
    pairs = result["matched_pairs"]
    constraint_best = front.sort_values(["constraint_score", "sequence_id"], ascending=[False, True]).iloc[0]
    novelty_best = front.sort_values(["novelty_score", "sequence_id"], ascending=[False, True]).iloc[0]
    conservative = front.loc[front.conservation_score.eq(1)].sort_values("constraint_score", ascending=False).iloc[0]
    feature_findings = f"""具体结果回答如下：

1. **前沿序列整体偏长，但前沿内仍有差异。**长度均值为 {primary.loc['length','front_mean']:.1f}，其余候选为 {primary.loc['length','control_mean']:.2f}；前沿范围为 {front.length.min():.0f}–{front.length.max():.0f} 个残基。长度匹配后的均值差降至 {matched.loc['length','mean_difference']:.1f} 个残基。
2. **组成更分散的描述在这一长度匹配对照中仍保留。**组成熵均值差由 {primary.loc['shannon_entropy_bits','mean_difference']:+.5f} 变为 {matched.loc['shannon_entropy_bits','mean_difference']:+.5f} bits；最高单一残基比例差由 {100*primary.loc['max_residue_fraction','mean_difference']:+.3f} 变为 {100*matched.loc['max_residue_fraction','mean_difference']:+.3f} 个百分点。它们描述整体组成分布，不能解释为功能更优。
3. **G/P 与同聚物差异依赖对照选择。**G/P 比例均值差由 {100*primary.loc['gly_pro_fraction','mean_difference']:+.3f} 变为 {100*matched.loc['gly_pro_fraction','mean_difference']:+.3f} 个百分点；最长同聚物段均值差由 {primary.loc['longest_homopolymer','mean_difference']:+.3f} 变为 {matched.loc['longest_homopolymer','mean_difference']:+.3f} 个残基，且两大组的中位数均为 {primary.loc['longest_homopolymer','front_median']:.0f}。不能把这两项表述为前沿共有的稳定优势。
4. **带电残基集合比例略低。**DEKR 比例均值差为 {100*primary.loc['charged_fraction','mean_difference']:+.3f} 个百分点，长度匹配后为 {100*matched.loc['charged_fraction','mean_difference']:+.3f} 个百分点。该组成描述不等于净电荷、可溶性或结构预测。

综上，当前前沿表现为不同目标取舍的序列集合；可观察到组成分散程度的组间差异，但没有证据把某一个组成特征称为全部前沿共享的功能机制。"""
    evidence = result["summary"]["evidence"]
    tie = result["tie_match_summary"].set_index("feature")
    tie_meta = result["summary"]["length_match_sensitivity"]
    tie_table = table(["特征", "100 解差值最小值", "中位数", "最大值", "负 / 零 / 正次数"], [
        [feature_names[f], f"{tie.loc[f,'min_mean_difference']:.5f}",
         f"{tie.loc[f,'median_mean_difference']:.5f}", f"{tie.loc[f,'max_mean_difference']:.5f}",
         f"{tie.loc[f,'negative_count']:.0f} / {tie.loc[f,'zero_count']:.0f} / {tie.loc[f,'positive_count']:.0f}"] for f in FEATURES])
    text = f'''# 选做实验 6：Pareto 前沿的目标取舍与序列特征

## 1. 实验问题与完成范围

本实验回答：四维非支配前沿中的候选分别在哪些维度突出；这些候选是否具有相似的序列特征。与实验 7 随机筛选对照共同构成两项选做，实验 5 下游性质预测未做。

当前固定排名池共有 114 条序列，第一前沿 10 条，其余 104 条。四维定义、资格门槛和原始评分均保持不变。本实验从发布清单核验后的分数与 FASTA 计算特征，不重新训练生成器。

## 2. 分析方案

分析前固定六项序列描述：长度、组成 Shannon 熵、最高单一残基比例、DEKR 比例、G/P 比例和最长同聚物段，并完整报告 20 种标准氨基酸组成。DEKR 为明确残基集合的计数比例，不是实际净电荷预测。

主比较为前沿 10 条与非前沿 104 条。敏感性对照使用最小总绝对长度差的一对一匹配，为 10 条前沿选取 10 条不同的非前沿候选。匹配目标仅包含长度，输入按 ID 排序；有多个等价最优解时使用固定软件环境得到的解，不把该解描述为唯一科学选择。实现使用 SciPy 的 [linear_sum_assignment](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linear_sum_assignment.html)。

组间给出均值、中位数、四分位数、范围及效应量。Cliff delta 为 P(前沿值大于对照) 减去 P(前沿值小于对照)。配对分析另保留每对的特征差值；这些特征没有统一的“越高越好”方向。所有比较为描述性分析，不计算 p 值，不进行功能或因果推断。

## 3. 前沿成员在哪些目标上突出？

![四维前沿](figures/frontier_objectives.png)

上图每条线是一条候选在四个原始目标上的得分，灰色为 104 条非前沿，绿色为四维第一前沿的 10 条；轴之间的连线只帮助阅读四维取舍。下图颜色表示各维度在同一 114 条候选池中的平均并列排名百分位；单元格为原始分数。百分位是相对位置，不是功能概率。最强/最弱维度按各维百分位比较，完整并列结果均保留。

{detail}

三个例子说明前沿不存在统一的最优序列：

- `{constraint_best.sequence_id}` 的约束分为 {constraint_best.constraint_score:.4f}，保守性为 {constraint_best.conservation_score:.4f}，体现较强的域约束支持。
- `{novelty_best.sequence_id}` 的新颖性为 {novelty_best.novelty_score:.4f}、独特性为 {novelty_best.uniqueness_score:.4f}，但约束分只有 {novelty_best.constraint_score:.4f}、保守性为 {novelty_best.conservation_score:.4f}，体现探索与约束之间的取舍。
- `{conservative.sequence_id}` 同时具有 {conservative.constraint_score:.4f} 的约束分和 {conservative.conservation_score:.4f} 的保守性；其新颖性为 {conservative.novelty_score:.4f}，低于上述探索端候选。

“最低百分位维度”只表示该候选各维相对排名中较低的一项，不表示不合格或生物学缺陷。例如保守性为 1 的候选仍可能因大量满分并列而在这列呈现较低百分位。

前沿按四个目标共同定义，不能把某个二维投影中的高点当成新的四维前沿。目标得分突出是前沿筛选规则的直接结果，不算独立验证。逐条可靠配对比例、预警与复核标记保存在 [front_details.tsv](front_details.tsv)。

## 4. 是否共享序列特征？

![六项特征](figures/sequence_features.png)

图中每点是一条候选，竖线为四分位距，横线为中位数；随机横向抖动只用于分开重叠点。匹配组是完整非前沿组的子集，不是第三组独立重复。以下均值差的方向统一为前沿减对照，完整分布摘要见 [group_comparison.tsv](group_comparison.tsv)。

{features_table}

{feature_findings}

![完整氨基酸组成](figures/amino_acid_composition.png)

{aa_table}

上述 20 种残基全部报告，不能只根据最大差异挑出一个残基并称其为决定功能的机制。只有 10 条前沿，且候选共享同一个生成来源；共有长度或组成范围可能来自训练分布及既有筛选门槛，不等于前沿特异的功能模式。

## 5. 长度匹配敏感性与证据覆盖

完成 {len(pairs)} 对不同候选的一对一匹配。绝对长度差总和为 {pairs.length_gap.sum():.0f} 个残基，中位数为 {pairs.length_gap.median():.1f}，最大为 {pairs.length_gap.max():.0f}。逐对对应关系见 [matched_pairs.tsv](matched_pairs.tsv)，逐对特征差见 [matched_feature_differences.tsv](matched_feature_differences.tsv)。匹配只降低长度差异，不消除其他组成、同源或缺失距离偏差。

### 等价最优匹配是否改变解释？

首次比较后发现，长度代价相同的最优解可能不唯一，因此补充固定 seed=42 的 100 次次级代价扰动。代价为“绝对长度差 × (前沿数+1) + [0,1) 随机数”，每次重新求全局一对一最优解，并断言原始总长度差仍为 {tie_meta['baseline_total_gap']}。主目标的一单位差大于全部次级成本之和，因而不会用更差的长度匹配换取某个特征方向。

得到 {tie_meta['unique_assignments']} 种配对、{tie_meta['unique_control_sets']} 种对照集合。该检查是事后补充的等价对照选择敏感性分析，不是 Bootstrap、显著性检验或 100 次独立生物学重复，也没有穷举全部最优解。以下差值全部为前沿减对照，比例特征保持 0–1 单位：

{tie_table}

组成熵差在这 100 个解中均为正，范围为 {tie.loc['shannon_entropy_bits','min_mean_difference']:.5f}–{tie.loc['shannon_entropy_bits','max_mean_difference']:.5f} bits；最高单一残基比例差均为负。G/P 差值 {tie.loc['gly_pro_fraction','negative_count']:.0f} 次为负、{tie.loc['gly_pro_fraction','positive_count']:.0f} 次为正，进一步说明它不是本次比较中方向稳定的前沿特征。原始主比较与匹配比较方向不同的结果全部保留，不能选择有利的一种解作结论。

全部匹配解见 [tie_match_pairs.tsv](tie_match_pairs.tsv)，全部特征差见 [tie_match_differences.tsv](tie_match_differences.tsv)。匹配不代表随机分组，也不证明长度是因果解释。

前沿仍有 {evidence['front']['manual_review_counts']['c_manual_review']}/10 条带人工复核标记。逐候选可靠配对比例均值为 {evidence['front']['uniqueness_resolved_fraction']['mean']:.2%}，其余候选为 {evidence['rest']['uniqueness_resolved_fraction']['mean']:.2%}，长度匹配对照为 {evidence['length_matched']['uniqueness_resolved_fraction']['mean']:.2%}。这里统计的是到原 136 条合格池的逐候选覆盖，不能混同入选子集内部的配对覆盖。前沿的探索性伴随更少的可解析距离证据，仍需保留复核。

完整摘要见 [summary.json](summary.json)；全部 13 层的描述见 [layer_summary.tsv](layer_summary.tsv)。距离未解析仍保留缺失状态，没有转成最大距离来奖励候选。

## 6. 复现与结论边界

```bash
python experiments/run_pareto_analysis.py --output results/reproductions/pareto_analysis
```

输出目录必须为新的或空目录。原 44 个评价结果不改变；新产物由独立清单记录。原始参考仍为暂定版本，当前只有一个生成批次，原生成权重未交付。本分析完成课程要求中的前沿目标解释与序列特征分析，不是下游性质预测或真实生物学功能验证。

源评价清单 SHA-256：`{source_hash}`。表和图由同一脚本生成，实际数值以 TSV 为准。
'''
    (output / "EXPERIMENT_6_REPORT.md").write_text(text, encoding="utf-8", newline="\n")


def build(source, output, qa_scripts=None):
    source, output = source.resolve(), output.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Output directory must be new or empty")
    manifest = json.loads((source / "manifest.json").read_text(encoding="utf-8"))
    verify_hashes(source, manifest["output_sha256"])
    scores = pd.read_csv(source / "ranking_pool_scores.tsv", sep="\t")
    sequences = read_fasta_strict(source / "ranking_pool.fasta")
    result = analyze_pareto_features(scores, sequences)
    ranking = pareto_ranking(scores, METRICS)
    expected = pd.read_csv(source / "pareto_top20.tsv", sep="\t")
    if ranking.head(20).sequence_id.tolist() != expected.sequence_id.tolist():
        raise ValueError("Recomputed four-objective ranking differs from published selection")
    if len(scores) != 114 or len(result["front_details"]) != 10:
        raise ValueError("This report describes the frozen 114-candidate pool and its 10-member front")
    output.mkdir(parents=True, exist_ok=True)
    for key, value in result.items():
        if isinstance(value, pd.DataFrame):
            value.to_csv(output / f"{key}.tsv", sep="\t", index=False, lineterminator="\n", float_format="%.12g")
    write_json(output / "summary.json", result["summary"])
    make_figures(result, output, qa_scripts)
    write_report(result, output, sha256_file(source / "manifest.json"))
    write_json(output / "manifest.json", {
        "source_manifest_sha256": sha256_file(source / "manifest.json"),
        "input_sha256": {f: sha256_file(source / f) for f in ["ranking_pool_scores.tsv", "ranking_pool.fasta", "pareto_top20.tsv"]},
        "implementation_sha256": {p: sha256_file(ROOT / p) for p in ["experiments/run_pareto_analysis.py", "src/gv_eval/pareto_features.py", "src/gv_eval/selection.py"]},
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "output_sha256": {p.relative_to(output).as_posix(): sha256_file(p) for p in sorted(output.rglob("*")) if p.is_file() and "qa" not in p.relative_to(output).parts},
    })
    print(json.dumps({"output": str(output), "front_count": len(result["front_details"]), "summary": result["summary"]}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--qa-script-dir", type=Path, help="Optional local figure-audit tools; not required to reproduce scientific data")
    args = parser.parse_args()
    build(args.source, args.output, args.qa_script_dir)
