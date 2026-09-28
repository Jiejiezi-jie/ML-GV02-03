# 复现与版本管理

## 环境

最小评价依赖固定在 `requirements-frozen-v2.txt`；QC、全局比对及证据校验使用 `requirements-integration.txt`。本机验证环境为 Python 3.12.14，VAE 测试使用 PyTorch 2.8.0 CPU。GPU 训练可在匹配 PyTorch CUDA 的独立环境执行。

仓库使用 LF 行尾。已经发布的原始证据目录在 `.gitattributes` 中保留精确字节，以免 Windows 检出改变哈希。

## 三种验证范围

| 命令 | 实际执行 | 输入条件 |
| --- | --- | --- |
| `experiments/run_full_experiment.py --output <新目录>` | 四维评分、九份名单、消融、鲁棒性、随机基线、两次输出一致性检查 | 最小依赖、完整 Git 历史 |
| `experiments/reproduce_project.py --output <新目录>` | 上述评价，另加两次 QC/全局比对、复核表汇总、与发布数据比较 | 扩展依赖、完整 Git 历史 |
| `experiments/check_frozen_checkout.py --report <新文件>` | 独立 Git 副本中的评价复现和全套测试 | 扩展依赖及 PyTorch |

独立副本检查使用当前选择的 Python 解释器及已安装依赖，清除 PYTHONPATH 和用户 site-packages，禁用外部 pytest 插件。它隔离源代码，不重新安装依赖。

同环境重复要求全部输出字节相同。跨环境比较先验证原清单，再核对 ID、行序、NaN 位置与数值（rtol=1e-9、atol=1e-12）；图片渲染的字节差异单独记录。

`--verify-existing` 是同代码、配置、环境下的严格验证，包含版本与图片哈希。不同环境核验发布批次使用扩展复现入口。

## 冻结来源

评价配置 `configs/d_evaluation.json` 固定生成批次与比对证据的源提交：

- 生成及初始评分：`10d9da8ef31cd8ec998ab53bec84252f331c64a6`。
- 全局比对及统一复核：`894f1fa80a8a0f5854e3816a9298c0c320197dbf`。

输入通过 `git show <commit>:<path>` 读取并校验哈希。所有源提交都已合入 main；分支名称是否保留不影响复现。下载 ZIP 或缺失这些对象的浅克隆不能运行冻结入口。修改工作区数据不会自动变更研究批次。

历史目录中的 `member_a_v1` 等名称属于已冻结路径标识。保留它们可避免改写来源清单和混淆批次，使用时以本仓库的正式入口及说明为准。

## 证据边界

扩展入口重算候选 QC、全局参考比对、同类距离、复核汇总和筛选实验。参考家族产物、全 Pfam 与 PureseqTM 发布表经过哈希核验，但该入口不重新调用这些外部工具。

原始 VAE checkpoint 应位于 `models/generator/sequence_vae/member_a_v1_seed42/best.pt`，当前未交付。训练和生成的记录保留原权重 SHA-256。重新训练会建立新权重、新批次，不能自动视作原始 1,000 条候选的重放。

集成复现原始记录位于 `results/gv02_03_v2/integration_20260928/`。其中 304 项测试及源码哈希属于清理前快照；交付整理后的检查见 `results/validation/`。

## 报告再生成

```bash
python experiments/build_project_report.py --output results/reproductions/report
```

入口验证冻结评价产物后，生成散点矩阵、相关图、策略图、鲁棒性图、Pareto 全池排名及 Markdown 报告。发布版本位于 `reports/`。演示材料内容同源于 `reports/statistics.json`。

`experiments/build_presentation.mjs` 保留演示的生成源码，使用 `@oai/artifact-tool` 和演示文稿技能运行时。运行时需设置 `ARTIFACT_TOOL_MODULE`、`RUNTIME_NODE_MODULES`、`PRESENTATIONS_SKILL_DIR`、`PRESENTATION_PYTHON` 为实际安装路径，再通过 Node.js 执行脚本。发布 PPTX 的图表和表格可直接在 PowerPoint 编辑；查看项目与复现实验不依赖此演示生成运行时。

## 历史材料

早期 200 条混合候选实验、旧演示、开发交接文档、原始大规模比对中间文件和旧生成模型已退出当前目录。完整清理前版本保存为 Git 标签 `archive/pre-delivery-cleanup-20260928`，对应提交 `e72b1da98fe1df6dce034ae8273feb5ec76475aa`。

```bash
git show archive/pre-delivery-cleanup-20260928:reports/M5_最终报告.md
```

原始数据和所有当前清单引用的模型、输入及输出继续保留。重建家族模型需要 Linux/WSL 中的 HMMER、MAFFT、CD-HIT 和 BLAST+，详见方法说明。
