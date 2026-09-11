# Fasta序列处理与分析项目

对应于T05任务中的输入fasta文件就是处理数据文件夹里面的deduplicated.fasta  # 去重后的fasta文件

## 项目概述

本项目包含一系列用于处理和分析fasta序列文件的Python脚本，主要功能包括：
1. 整合多个文件夹中的fasta文件
2. 基于90%相似性阈值对序列进行去重
3. 对去重后的序列进行统计分析，包括序列长度分布和氨基酸频率

## 目录结构

```
├── 原始数据/               # 原始fasta文件存放目录
│   ├── jtl/               # 子目录1
│   ├── jxzw/              # 子目录2
│   ├── ljl/               # 子目录3
│   ├── lxy/               # 子目录4
│   └── sy/                # 子目录5
├── 汇总数据/               # 整合后的fasta文件和脚本
│   ├── merged.fasta        # 整合所有fasta文件后的结果
│   └── merge_fasta.py      # 整合fasta文件脚本
├── 处理数据/               # 去重后的fasta文件和脚本
│   ├── deduplicated.fasta  # 去重后的fasta文件
│   └── deduplicate_fasta.py # 序列去重脚本
├── 分析数据/               # 统计分析结果和脚本
│   ├── analyze_fasta.py    # 序列统计分析脚本
│   ├── amino_acid_frequency.png  # 氨基酸频率柱状图
│   ├── amino_acid_frequency.txt  # 氨基酸频率统计
│   ├── length_boxplot.png        # 序列长度箱线图
│   ├── length_distribution.png   # 序列长度分布直方图
│   ├── length_distribution.txt   # 所有序列长度列表
│   └── length_stats.txt          # 序列长度统计指标
└── README.md              # 项目说明文档
```

## 脚本功能说明

### 1. merge_fasta.py

**功能**：整合原始数据目录及其子目录下所有的fasta文件到一个文件中。

**使用方法**：
```bash
python 汇总数据/merge_fasta.py
```

**生成文件**：
- `汇总数据/merged.fasta`：整合了所有fasta文件的序列

### 2. deduplicate_fasta.py

**功能**：基于CD-HIT算法思路，对fasta序列进行去重，保留90%相似性阈值下的代表性序列。

**使用方法**：
```bash
python 处理数据/deduplicate_fasta.py
```

**参数说明**：
- 输入文件：默认`汇总数据/merged.fasta`
- 输出文件：默认`处理数据/deduplicated.fasta`
- 相似性阈值：默认90%

**生成文件**：
- `处理数据/deduplicated.fasta`：去重后的序列文件

### 3. analyze_fasta.py

**功能**：对fasta序列进行统计分析，包括：
- 序列长度分布（最小值、最大值、平均值、中位数、标准差）
- 氨基酸频率统计（作为"真实分布"基准）
- 生成可视化图表

**使用方法**：
```bash
python 分析数据/analyze_fasta.py
```

**参数说明**：
- 输入文件：默认`处理数据/deduplicated.fasta`

**生成文件**：
- 文本文件：
  - `分析数据/length_distribution.txt`：所有序列长度列表
  - `分析数据/length_stats.txt`：序列长度统计指标
  - `分析数据/amino_acid_frequency.txt`：氨基酸频率统计（真实分布基准）
- 可视化图表：
  - `分析数据/length_distribution.png`：序列长度分布直方图
  - `分析数据/amino_acid_frequency.png`：氨基酸频率分布柱状图
  - `分析数据/length_boxplot.png`：序列长度箱线图

## 统计结果说明

### 序列长度统计

| 统计指标 | 数值 |
|----------|------|
| 序列数量 | 4379条 |
| 最短序列 | 54 aa |
| 最长序列 | 757 aa |
| 平均长度 | 148.02 aa |
| 中位数长度 | 106.0 aa |
| 长度标准差 | 95.28 aa |

### 氨基酸频率（真实分布基准，前10位）

| 氨基酸 | 频率 | 出现次数 |
|--------|------|----------|
| L      | 12.81% | 83,049次 |
| E      | 9.86%  | 63,921次 |
| R      | 8.41%  | 54,534次 |
| D      | 8.19%  | 53,110次 |
| A      | 7.93%  | 51,382次 |
| V      | 7.59%  | 49,185次 |
| G      | 6.86%  | 44,444次 |
| T      | 5.71%  | 37,026次 |
| S      | 5.31%  | 34,428次 |
| I      | 4.85%  | 31,426次 |

## 依赖项

- Python 3.x
- Biopython (`pip install biopython`)
- NumPy (`pip install numpy`)
- Matplotlib (`pip install matplotlib`)

## 使用流程

1. **整合fasta文件**：运行`汇总数据/merge_fasta.py`，生成`汇总数据/merged.fasta`
2. **序列去重**：运行`处理数据/deduplicate_fasta.py`，生成`处理数据/deduplicated.fasta`
3. **统计分析**：运行`分析数据/analyze_fasta.py`，生成统计结果和图表到`分析数据/`目录

## 注意事项

1. 脚本分布在不同的子目录中，运行时需使用正确的相对路径
2. 脚本会自动处理路径关系，无需手动调整输入输出路径
3. 确保原始数据目录下有fasta文件，否则整合脚本将无法生成数据
4. 生成的图表可能会显示中文字符警告，但不影响图表内容
5. 建议使用UTF-8编码打开生成的文本文件

## 结果解释

- **序列长度分布**：反映了序列的长度特征，有助于了解蛋白质的大小范围
- **氨基酸频率**：作为"真实分布"基准，可用于后续的序列比对、进化分析等
- **可视化图表**：直观展示了序列特征，便于快速理解数据分布

## 扩展功能

可根据需要修改脚本中的参数：
- 调整去重的相似性阈值
- 修改图表的样式和大小
- 添加更多统计分析指标

