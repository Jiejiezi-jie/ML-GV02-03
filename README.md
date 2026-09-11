# GV02-03：GvpA 候选序列的多目标评价与筛选

本项目属于教师指南 [project-pool(1).pdf](project-pool(1).pdf) 的 **Track 1：GV（气囊蛋白）／GV02-03**。我们接收已有候选，评价域约束满足度、关键位点保守性、新颖性和多样性，并比较筛选策略、开展消融及鲁棒性分析。多样性是四个目标之一；本项目不要求训练新生成模型。

**候选入口：[T05 生成候选集](data/candidates/t05/generated_gvp.fasta)，共 200 条。** 这些序列来自混合 Gvp 数据训练的生成模型，尚未被核验为纯 GvpA。天然序列用于建立参考集，不能当作生成候选混入评价。

## 当前阶段

本轮完成 M1 问题分析：理解背景、分析数据特点、抽象和定义任务，形成 [功能需求分析与任务建模报告](reports/M1_功能需求分析与任务建模报告.md)（[PDF](reports/M1_功能需求分析与任务建模报告.pdf)）。输入统计是实际审计结果；域扫描、序列比对、正式四目标评分和筛选实验属于后续工作。

[机器学习要求.md](机器学习要求.md) 保留原始任务记录。任务范围以课程指南的 GV02-03 为准，T05 原报告用于追溯输入来源。

## 项目目录

```text
machine/
├── project-pool(1).pdf              # 教师课程指南
├── 机器学习要求.md                  # 原始任务记录
├── data/
│   ├── README.md                   # 六份 FASTA 的用途和风险
│   ├── candidates/t05/             # 200 条生成候选
│   ├── raw/t05/gvpa/               # 三份名义 GvpA 天然来源
│   ├── raw/t05/real_gvp.fasta       # 混合 Gvp 来源和训练背景
│   └── raw/design/GvpA.fasta        # 有独有序列的天然参考来源
├── references/t05/                # 原始报告、方法说明及生成资源
├── preprocessing/                 # 输入审计脚本
├── results/input_audit/            # 当前统计摘要和逐候选审计表
├── reports/                        # M1 报告
└── docs/
    ├── 整理说明.md                  # 目录整理说明
    └── cleanup/                    # 两轮清理与校验记录
```

第二轮删除 32 个无关或重复文件，共 6,811,715 字节。design 工程、design 的 10 条旧候选和评分、T05 重复演示材料及旧处理/分析脚本已删除。design 天然 FASTA 因有 30 条其余三份天然来源未覆盖的完整序列而保留。详见 [整理说明](docs/整理说明.md)。

## 数据风险

- 当前有 **6 个 FASTA 文件**。四份名义 GvpA 来源的不同完整序列并集为 1252 条，尚未完成家族和质量过滤。
- T05 候选长 66–512 aa，中位数 103 aa，全部使用标准氨基酸且无完整序列重复。两条长度恰为生成上限 512 aa，需要核查是否截断。
- 天然来源存在重复、未知残基 `X`、不完整注释和 GvpJ 等其他家族条目。不能仅按长度判断 GvpA 身份；没有精确重复也不能证明高新颖性。
- 原任务记录中的 **PF01132 是 EFP（EF-P 的 OB 域），不能作为 GvpA 域标识**。相关的 PF00741／Gas_vesicle 家族也包含 GvpJ，因此命中不足以单独确定 GvpA 亚家族；后续需结合参考注释、比对与覆盖度核查，并固定数据库版本。[NCBI PF01132](https://www.ncbi.nlm.nih.gov/Structure/cdd/pfam01132)、[NCBI PF00741](https://www.ncbi.nlm.nih.gov/Structure/cdd/pfam00741)

## 复核输入

在项目根目录运行（Python 3.8 及以上，无第三方依赖）：

```powershell
python -B -X utf8 preprocessing/analyze_inputs.py
```

命令读取原始 FASTA，更新 [summary.json](results/input_audit/summary.json) 和 [candidate_records.tsv](results/input_audit/candidate_records.tsv)，不修改原始序列，不执行域扫描、比对或家族判定。具体来源见 [数据说明](data/README.md)，原生成资源的缺失依赖见 [参考材料说明](references/README.md)。

`docs/cleanup/` 第一轮 `fasta_inventory.json` 等文件是历史快照，包含当时的七份 FASTA；**当前统计以 `results/input_audit/summary.json` 为准**。后续清洗参考集、评分表、筛选名单和实验结果应另行输出，保持原始输入不变。

报告正文可直接编辑 Markdown。重新生成 PDF 和长度分布图需要 Python 的 `reportlab`、Poppler 的 `pdftoppm`，以及 Windows 宋体和黑体字体：运行 `python -B reports/build_report.py --pdftoppm "pdftoppm.exe的实际路径"`。输入审计无需这些额外依赖。
