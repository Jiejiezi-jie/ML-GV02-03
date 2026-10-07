# 复现与版本管理

## 环境

最小评价及前沿特征分析依赖固定在 `requirements-frozen-v2.txt`；QC、全局比对及证据校验使用 `requirements-integration.txt`。2026-09-28 发布验证环境为 Python 3.12.14，VAE 测试使用 PyTorch 2.8.0 CPU；各次运行的实际版本以相应清单为准。GPU 训练可在匹配 PyTorch CUDA 的独立环境执行。

仓库使用 LF 行尾。已经发布的原始证据目录在 `.gitattributes` 中保留精确字节，以免 Windows 检出改变哈希。

## 验证范围

| 命令 | 实际执行 | 输入条件 |
| --- | --- | --- |
| `experiments/run_full_experiment.py --output <新目录>` | 四维评分、九份名单、消融、鲁棒性、随机基线、两次输出一致性检查 | 最小依赖、完整 Git 历史 |
| `experiments/run_pareto_analysis.py --output <新目录>` | 核验发布评价后，计算四维前沿、序列特征、长度匹配、等价匹配敏感性及专项报告 | 最小依赖、完整的发布评价结果 |
| `experiments/reproduce_project.py --output <新目录>` | 上述评价，另加两次 QC/全局比对、复核表汇总、与发布数据比较 | 扩展依赖、完整 Git 历史 |
| `experiments/check_frozen_checkout.py --report <新文件>` | 独立 Git 副本中的评价复现和全套测试 | 扩展依赖及 PyTorch |

独立副本检查使用当前选择的 Python 解释器及已安装依赖，清除 PYTHONPATH 和用户 site-packages，禁用外部 pytest 插件。它隔离源代码，不重新安装依赖。

同环境重复要求全部输出字节相同。跨环境比较先验证原清单，再核对 ID、行序、NaN 位置与数值（rtol=1e-9、atol=1e-12）；图片渲染的字节差异单独记录。

`--verify-existing` 是同代码、配置、环境下的严格验证，包含版本与图片哈希。不同环境核验发布批次使用扩展复现入口。

## Pareto 前沿分析

```bash
python experiments/run_pareto_analysis.py --output results/reproductions/pareto_analysis
```

该入口先核验 `results/gv02_03_v2/d_evaluation_member_a_v1/manifest.json` 中的 44 个发布输出，再从 `ranking_pool_scores.tsv`、`ranking_pool.fasta` 重算前沿并核对发布 Top-20。分析保持 114 条排名池和原始四维分数不变；输出目录必须为新的或空目录。

选做实验 6 的发布目录为 `results/gv02_03_v2/pareto_analysis_member_a_v1/`，包含：

- 六项预定特征和全部 20 种氨基酸频率、10 条前沿逐候选解释及全部 13 层摘要。
- 10 条前沿对 104 条非前沿的完整描述；最小总长度差的一对一匹配及逐对差值。
- 固定种子 42 的 100 次等价最优匹配检查：原始总长度差固定为 2，以次级随机代价检查多个最优解的影响。
- Markdown 报告、图表，以及记录来源、代码、依赖版本和输出哈希的独立 `manifest.json`。

一次调用生成一次结果；它不自动执行第二次比较。核对同环境可重复性时，使用两个空目录分别运行，再按清单逐项比较输出哈希。图表受软件及字体版本影响，跨环境需先区分数值结果与渲染字节。匹配敏感性不是 Bootstrap 或独立实验重复，全部特征比较保持描述性解释。

2026-10-07 在同一环境中独立运行两次，24 个表格、报告和图表产物及其清单全部逐字节一致；核验记录见 [experiment6_20261007](../results/validation/experiment6_20261007/README.md)。本次全套 333 项测试通过。Windows 运行测试时设置 `PYTHONUTF8=1`、`PYTHONIOENCODING=utf-8`，使子进程输出也使用 UTF-8，避免系统默认编码导致日志解码警告。

三份新图同时提供可编辑文本的 SVG/PDF 和 300 dpi PNG。图表检查记录位于实验目录的 `qa/`，检查工具不属于运行实验的必需依赖。

## 冻结来源

评价配置 `configs/d_evaluation.json` 固定生成批次与比对证据的源提交：

- 生成及初始评分：`10d9da8ef31cd8ec998ab53bec84252f331c64a6`。
- 全局比对及统一复核：`894f1fa80a8a0f5854e3816a9298c0c320197dbf`。

输入通过 `git show <commit>:<path>` 读取并校验哈希。所有源提交都已合入 main；分支名称是否保留不影响复现。下载 ZIP 或缺失这些对象的浅克隆不能运行冻结入口。修改工作区数据不会自动变更研究批次。

历史目录中的 `member_a_v1` 等名称属于已冻结路径标识。保留它们可避免改写来源清单和混淆批次，使用时以本仓库的正式入口及说明为准。

## 证据边界

扩展入口重算候选 QC、全局参考比对、同类距离、复核汇总和筛选实验。参考家族产物、全 Pfam 与 PureseqTM 发布表经过哈希核验，但该入口不重新调用这些外部工具。

原始 VAE checkpoint 应位于 `models/generator/sequence_vae/member_a_v1_seed42/best.pt`，当前未交付。训练和生成的记录保留原权重 SHA-256。重新训练会建立新权重、新批次，不能自动视作原始 1,000 条候选的重放。

清理前的 304 项测试及集成记录已归档，可在标签 `archive/pre-final-delivery-20261007` 下的 `results/gv02_03_v2/integration_20260928/` 查看。2026-09-28 交付整理后的 302 项测试、完整复现与独立副本检查保留在 `results/validation/`；2026-10-07 的前沿分析验证见其 `experiment6_20261007/` 子目录，333 项测试通过。原 44 个评价输出保持不变。历史记录保留原始日期和范围，不代表本次重新执行了外部建模或扫描。

## 报告再生成

```bash
python experiments/build_project_report.py --output results/reproductions/report
```

入口验证冻结评价及实验 6 的产物后，生成散点矩阵、相关图、策略图、鲁棒性图、Pareto 全池排名及包含实验 6 结论的 Markdown 报告。发布版本位于 `reports/`。实验 6 的专项报告与图表由 `run_pareto_analysis.py` 单独生成，正式报告已纳入该扩展。中期演示保留当时范围，结题答辩 PPT 尚未制作。

更新已发布报告使用 `python experiments/build_project_report.py --replace-published`。入口先校验已有清单，再仅更新八项生成报告产物及清单，保留演示稿等其他文件。正式报告清单不覆盖中期材料；重新生成报告不会改动中期演示。

### 中期演示与结题状态

中期答辩 PPT、逐页讲稿和已有项目概览保存在 [reports/midterm/](../reports/midterm/README.md)，原始文件字节保持不变。用户修改过的 13 页概览来自提交 `7482dcac254c30f84f2ea147300b3200c1ed8a84`，当前文件名为 `PROJECT_OVERVIEW.pptx`。

结题答辩 PPT 尚未制作。中期材料采用当时的实验范围，不作为涵盖选做 6 新分析的结题演示。历史演示生成源码可从 `archive/pre-final-delivery-20261007:experiments/build_presentation.mjs` 查看；数值实验和报告复现不依赖该演示运行时。

## 历史材料

早期 200 条混合候选实验、旧演示、开发交接文档、原始大规模比对中间文件和旧生成模型已退出当前目录。完整清理前版本保存为 Git 标签 `archive/pre-delivery-cleanup-20260928`，对应提交 `e72b1da98fe1df6dce034ae8273feb5ec76475aa`。

```bash
git show archive/pre-delivery-cleanup-20260928:reports/M5_最终报告.md
```

原始数据和所有当前清单引用的模型、输入及输出继续保留。重建家族模型需要 Linux/WSL 中的 HMMER、MAFFT、CD-HIT 和 BLAST+，详见方法说明。

课程评分与工程要求原件 `docs/course/course-introduction.pdf` 按原始字节从上述标签的 `intro-mlproj.pdf` 恢复；项目题目原件与过程要求分别位于同目录的 `project-pool.pdf`、`process-requirements.docx`。

2026-10-07 的目录整理另保存标签 `archive/pre-final-delivery-20261007`，对应 `8aa3c81`。清理前集成快照已退出当前目录，可通过该标签恢复。中期演示与讲稿保留在 `reports/midterm/`，不参与数值实验运行。已有 13 页项目概览的本地修改稿另保存在提交 `7482dcac254c30f84f2ea147300b3200c1ed8a84`。

```bash
git ls-tree -r --name-only archive/pre-final-delivery-20261007 reports/midterm
git show archive/pre-final-delivery-20261007:results/gv02_03_v2/integration_20260928/reproduction_report.json
```

当前仍保留的家族、域、比对和生成证据均被配置、清单或测试引用。目录名中的历史批次标识不影响其作为正式输入证据的用途，不能按名称批量删除。
