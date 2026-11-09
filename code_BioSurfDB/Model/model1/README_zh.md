# BioSurfDB 分子指纹建模 —— 任务一：通用生物学功能注释状态分类

## 项目简介

本任务使用 Morgan 分子指纹（radius=2，nBits=1024）与机器学习模型（随机森林、XGBoost），对化合物的**通用生物学功能注释状态**做二分类预测（`biological_activity`，YES / NO）。

- **阳性样本（pos）**：`biological_pos.xlsx`，其 `biological_activity` 列值为 `YES`。
- **比较样本（comp）**：`biological_comp.xlsx`，其 `biological_activity` 列值为 `NO`，作为对照的第二类。

> 说明：`biological_comp` 是"比较 / 参考"样本集，并非绝对"无生物活性"，仅作为二分类中的对照群体。驱动模型的标签始终来自 `biological_activity` 列。

采用 5 折分层交叉验证评估模型，输出 ROC 曲线、各项分类指标（Accuracy、Precision、Recall、F1、ROC-AUC、MCC、Balanced Accuracy）以及 Morgan 指纹位点重要性分析。

## 目录结构

```text
model1/
├── input/
│   ├── biological_pos.xlsx             # 阳性样本（生物学功能注释状态 = YES）
│   ├── biological_comp.xlsx            # 比较/参考样本（生物学功能注释状态 = NO）
│   ├── fingerprints_biological/        # 脚本 01 生成的分子指纹
│   │   ├── fingerprints_pos.csv        # 阳性样本指纹
│   │   └── fingerprints_comp.csv       # 比较样本指纹
│   └── cv_folds_biological/            # 脚本 02 生成的 5 折划分
│       ├── fold_1_train.csv
│       ├── fold_1_val.csv
│       └── ...
├── output/
│   ├── biological/                     # 任务一输出
│   │   ├── figures/
│   │   │   ├── rf_roc_curves.png
│   │   │   └── xgb_roc_curves.png
│   │   ├── rf_cv_results.csv
│   │   ├── rf_summary.csv
│   │   ├── xgb_cv_results.csv
│   │   └── xgb_summary.csv
│   └── bit_importance_biological/      # 位点重要性输出
│       ├── rf_stable_bits.png
│       ├── xgb_stable_bits.png
│       └── bit_importance_summary.csv
├── scripts/
│   ├── 01_generate_fingerprints.py     # 生成分子指纹
│   ├── 02_cross_validation_5fold.py    # 5 折分层数据划分
│   ├── 03_train_random_forest.py       # 随机森林训练
│   ├── 04_train_xgboost.py             # XGBoost 训练
│   └── 05_bit_importance_analysis.py   # 指纹位点重要性分析
├── environment.yml                     # Conda 最小核心依赖（推荐，搭建快）
├── environment_full.yml                # 完整锁定快照（锁定全部传递依赖，用于逐位复刻）
├── requirements.txt                    # pip 依赖（不含 RDKit）
└── README_en.md / README_zh.md
```

## 环境配置

### 方式一：使用 Conda 环境（推荐）

```bash
conda env create -f environment.yml
conda activate surfmodel
```

> **完全复刻：** 若需重建论文发表时使用的逐位一致环境（锁定全部传递依赖与构建字符串），请改用 `environment_full.yml` 创建环境：
> ```bash
> conda env create -f environment_full.yml
> conda activate surfmodel
> ```
> `environment_full.yml` 为 `conda env export` 快照，与平台相关（抓取于 Windows）；仅在需要严格可复现时启用。

- `environment.yml` 中的频道已**写全为清华镜像 URL**（`conda-forge` + `main`），自包含、可直接复现，不依赖任何全局 `.condarc` 配置，国内可直接使用。
- 若 pip 安装阶段（如 `xgboost`）偏慢，可先设置 PyPI 镜像：
  ```bash
  pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple
  ```

### 方式二：手动安装（备选）

```bash
pip install -r requirements.txt
conda install -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ rdkit
```

> 注意：RDKit 必须通过 conda 安装；`requirements.txt` 中已列出除 RDKit 外的所有包。

## 数据准备

`biological_pos.xlsx` 与 `biological_comp.xlsx` 已放在 `input/` 文件夹下，无需另行准备。

输入 Excel 需包含以下由上游 RDKit 脚本（`../RDkit/deal_RDkit.py`）标准化后生成的列，脚本实际读取其中 4 列：

| 列名 | 说明 |
| --- | --- |
| `Database_ID` | 化合物编号 |
| `SMILES_Parent_RDKit` | 脱盐/去电荷后的规范 SMILES，用于生成分子指纹（唯一指纹来源列） |
| `biological_activity` | **任务一标签**（YES/NO → 1/0） |
| `antimicrobial_activity` | **任务二标签**（抗菌/抑菌活性，当前任务一不使用，保留备用） |

> 标签映射规则（`map_label`）：`YES`/`TRUE`/`1` → 1，其余（含 `NO`、空值、`NaN`）→ 0。

## 运行步骤

所有脚本均通过 `Path(__file__).resolve().parent.parent` 定位项目根目录并使用绝对路径，可在任意位置运行。

```bash
conda activate surfmodel
cd Model/model1/scripts
```

1. 生成分子指纹

   ```bash
   python 01_generate_fingerprints.py
   ```

   输出：`input/fingerprints_biological/fingerprints_pos.csv`、`fingerprints_comp.csv`

2. 5 折分层交叉验证划分

   ```bash
   python 02_cross_validation_5fold.py
   ```

   输出：`input/cv_folds_biological/fold_*_train.csv`、`fold_*_val.csv`

3. 训练随机森林

   ```bash
   python 03_train_random_forest.py
   ```

   输出：`output/biological/rf_cv_results.csv`、`rf_summary.csv`、`figures/rf_roc_curves.png`

4. 训练 XGBoost

   ```bash
   python 04_train_xgboost.py
   ```

   输出：`output/biological/xgb_cv_results.csv`、`xgb_summary.csv`、`figures/xgb_roc_curves.png`

5. 指纹位点重要性分析

   ```bash
   python 05_bit_importance_analysis.py
   ```

   输出：`output/bit_importance_biological/` 下的图片与 `bit_importance_summary.csv`

> 提示：若之前用旧命名跑过脚本，请删除 `input/fingerprints_biological/` 下残留的 `fingerprints_positive.csv` / `fingerprints_negative.csv`，以免与新生成的 `fingerprints_pos.csv` / `fingerprints_comp.csv` 混淆。

## 结果解读

- **CSV 结果文件**：每个折叠的 Accuracy、Precision、Recall、F1、ROC-AUC、MCC、Balanced Accuracy。
- **汇总文件**（`*_summary.csv`）：各指标在 5 折上的均值。
- **ROC 曲线图**：曲线越靠近左上角，性能越好。
- **位点重要性图**：展示在至少 3 折中稳定出现的高重要性 Morgan 指纹位点（基于排列重要性）。

## 任务二（可选）

脚本 `05_bit_importance_analysis.py` 同时支持抗菌活性任务：将文件顶部的 `TASK` 改为 `"antimicrobial"`，标签列自动切换为 `antimicrobial_activity`，并读取 `input/cv_folds_antimicrobial/`（需另行准备对应的 5 折划分）。当前主线流程以任务一（`biological`）为准。

## 注意事项

- 所有脚本均通过 `Path(__file__).resolve().parent.parent` 定位项目根目录，使用绝对路径，无需手动 `cd` 到特定目录。
- 随机种子统一为 `RANDOM_STATE = 42`，结果可复现。
- 分子指纹来源列固定为 `SMILES_Parent_RDKit`。
- 指纹半径 `RADIUS = 2`、维度 `NBITS = 1024`；随机森林使用 `class_weight="balanced"` 处理类别不平衡。

## 常见问题

**Q：缺少 openpyxl？**
A：已包含在 `environment.yml` 和 `requirements.txt` 中，重新创建/安装环境即可。

**Q：RDKit 安装失败？**
A：使用 conda 单独安装（已配清华镜像）：

```bash
conda install -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/ rdkit
```

**Q：运行 01 时报 `KeyError` 或找不到列？**
A：检查输入 Excel 是否包含 `biological_activity` 与 `antimicrobial_activity` 列（名称需完全一致，区分大小写）。
