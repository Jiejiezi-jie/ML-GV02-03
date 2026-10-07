# 最终目录验收

日期：2026-10-07。验证对象为清理后的 GV02-03 项目。

## 实验与代码

- 使用独立 Git 副本运行 `experiments/check_frozen_checkout.py`，包括完整评价复跑及全套测试：333 项通过。
- 与原发布评价相比，44 个产物中的 42 个字节一致；`metric_distributions.png` 和 `strategy_diversity.png` 只有渲染字节差异。科学数据、行序和候选名单一致，全部比较通过。
- 在独立输出目录运行前沿分析：24 个产物和其清单均与发布文件逐字节一致。
- 保留两套正式实验清单，原始结果未被清理或改写。

原始独立副本记录见 `checkout_acceptance.json`。该运行在目录移除后、教学文档整理期间执行，因此其中 `snapshot_sha256` 是当时的文档快照；科学实现及输入未随后修改。最终文档、中期材料和链接另由 `delivery_validation.json` 核验。

## 目录与材料

旧集成快照已在 `archive/pre-final-delivery-20261007` 标签保留，当前工作树移除其副本；中期材料保留原版本并集中到 `reports/midterm/`。用户修改过的演示原稿还保存在 `7482dcac254c30f84f2ea147300b3200c1ed8a84`。最终交付不包含重复复跑目录或 Python/pytest 缓存。

`delivery_validation.json` 记录实际清理项目、最终材料哈希及核验结果；`manifest.json` 对本目录文件提供校验。结题答辩 PPT 尚未制作；真实过程参与和现场答辩不在自动验收范围内。
