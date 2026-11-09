# BioSurf SMILES 标准化工具（RDKit）

使用 RDKit 对化合物 SMILES 进行标准化处理，自动从 `SMILES_PubChem` 与 `SMILES_Original` 两列中选取 SMILES（优先 `SMILES_PubChem`，为空时回退 `SMILES_Original`），并生成规范化的化学信息列。

## 功能

对每一行 SMILES 执行：脱盐（保留母核片段）→ 去电荷（得到父结构）→ 生成规范 SMILES、InChI、InChIKey、分子式、精确分子量（单同位素）。两列均为空的行会跳过并计数。

输出新增列：

| 列名 | 含义 |
|---|---|
| `SMILES_RDKit` | 脱盐后的规范 SMILES |
| `SMILES_Parent_RDKit` | 脱盐并去电荷后的父结构规范 SMILES |
| `InChI_RDKit` | InChI 字符串 |
| `InChIKey_RDKit` | InChIKey |
| `Molecular_Formula_RDKit` | 分子式 |
| `Molecular_Weight_RDKit` | 精确分子量（单同位素，非平均分子量） |

## 目录结构

```
RDkit/
├── environment.yml        # Conda 环境定义（已配置清华镜像）
├── environment_full.yml    # 完整锁定快照（锁定全部传递依赖，用于逐位复刻）
├── requirements.txt       # 纯 pip 依赖（与 environment.yml 的 pip 段一致）
├── input/                 # 放置输入 Excel（运行前需自行放入）
├── output/                # 运行后自动生成结果（脚本会自动创建）
├── scripts/
│   └── deal_RDkit.py      # 主处理脚本
└── README_zh.md / README_en.md
```

## 环境搭建

需要已安装 Anaconda 或 Miniconda。

```bash
# 在项目根目录（RDkit/）执行
conda env create -f environment.yml
conda activate surfdb
```

> **完全复刻：** 若需重建论文发表时使用的逐位一致环境，请改用 `environment_full.yml` 创建环境：
> ```bash
> conda env create -f environment_full.yml
> conda activate surfdb
> ```
> `environment_full.yml` 为 `conda env export` 快照，与平台相关（抓取于 Windows）；仅在需要严格可复现时启用。

- 若 pip 安装 `openpyxl` / `pubchempy` 较慢，先设置 PyPI 镜像：
  ```bash
  pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple
  ```

> 已有 conda 环境、只想补装纯 Python 依赖时，也可：`pip install -r requirements.txt`
> （注意 `rdkit` 无法通过 pip 可靠安装，仍需用 conda 安装）。

## 使用方法

1. 输入 Excel 已放在 `input/` 目录中（默认文件名 `data_table_Biologic_core.xlsx`）。
   - 若文件名不同，修改 `scripts/deal_RDkit.py` 中第 174 行左右的 `input_file` 文件名即可。
   - 输入表需包含 `SMILES_PubChem` 或 `SMILES_Original` 列（按子串匹配，脚本读取第一张工作表）。
2. 运行（任意目录下均可，脚本按自身位置自动定位 input/output）：
   ```bash
   python scripts/deal_RDkit.py
   ```
3. 结果输出到 `output/data_table_Biologic_core_standardized.xlsx`。

运行结束后会打印统计：总化合物数、成功处理数、失败数、两列均为空的行数。

## 依赖版本

- python 3.10.20
- rdkit 2026.03.2
- pandas 2.3.3
- numpy 2.2.6
- openpyxl 3.1.5（pip）
- pubchempy 1.0.5（pip）
