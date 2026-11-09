# BioSurf SMILES Standardisation Tool (RDKit)

Standardises compound SMILES using RDKit. It automatically picks a SMILES from the
`SMILES_PubChem` and `SMILES_Original` columns (preferring `SMILES_PubChem`, falling
back to `SMILES_Original` when empty) and produces normalised chemical-information columns.

## Features

For each SMILES row it performs: desalting (keep the parent fragment) → charge
neutralisation (parent structure) → generate canonical SMILES, InChI, InChIKey,
molecular formula, and exact molecular weight (monoisotopic). Rows where both SMILES
columns are empty are skipped and counted.

New output columns:

| Column | Meaning |
|---|---|
| `SMILES_RDKit` | Canonical SMILES after desalting |
| `SMILES_Parent_RDKit` | Canonical SMILES of the desalted & uncharged parent |
| `InChI_RDKit` | InChI string |
| `InChIKey_RDKit` | InChIKey |
| `Molecular_Formula_RDKit` | Molecular formula |
| `Molecular_Weight_RDKit` | Exact molecular weight (monoisotopic, not average) |

## Directory layout

```
RDkit/
├── environment.yml        # Conda environment definition (Tsinghua mirror configured)
├── environment_full.yml    # Full locked snapshot (all transitive deps pinned; exact reproduction)
├── requirements.txt       # Pure pip dependencies (matches the pip section of environment.yml)
├── input/                 # Place the input Excel here (add it before running)
├── output/                # Results are written here (auto-created by the script)
├── scripts/
│   └── deal_RDkit.py      # Main processing script
└── README_zh.md / README_en.md
```

## Environment setup

Requires Anaconda or Miniconda installed.

```bash
# Run from the project root (RDkit/)
conda env create -f environment.yml
conda activate surfdb
```

> **Exact reproduction:** to recreate the bit-identical environment used at publication time, create it from `environment_full.yml` instead:
> ```bash
> conda env create -f environment_full.yml
> conda activate surfdb
> ```
> `environment_full.yml` is a `conda env export` snapshot and is platform-specific (captured on Windows); use it only when strict reproducibility is required.

- If `openpyxl` / `pubchempy` install slowly via pip, set a PyPI mirror first:
  ```bash
  pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple
  ```

> To add only the pure-Python deps to an existing conda env: `pip install -r requirements.txt`
> (note: `rdkit` cannot be reliably installed via pip and must come from conda).

## Usage

1. The input Excel is already provided in `input/` (default name `data_table_Biologic_core.xlsx`).
   - For a different filename, edit the `input_file` name near line 174 of `scripts/deal_RDkit.py`.
   - The input sheet must contain a `SMILES_PubChem` or `SMILES_Original` column
     (matched by substring; the first worksheet is read).
2. Run (from any directory — the script locates input/output relative to itself):
   ```bash
   python scripts/deal_RDkit.py
   ```
3. Results are written to `output/data_table_Biologic_core_standardized.xlsx`.

After finishing, it prints statistics: total compounds, successfully processed,
failed, and rows where both SMILES columns were empty.

## Dependency versions

- python 3.10.20
- rdkit 2026.03.2
- pandas 2.3.3
- numpy 2.2.6
- openpyxl 3.1.5 (pip)
- pubchempy 1.0.5 (pip)
