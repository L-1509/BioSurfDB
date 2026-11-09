# BioSurfDB — Code Repository

This repository contains the code used to build the **BioSurfDB** core dataset and the
machine-learning analysis presented in the accompanying paper.

## What's inside

| Module | What it does |
|---|---|
| `code_BioSurfDB/RDkit/` | Chemical-structure standardisation with RDKit. Reads the raw compound tables, desalinates / neutralises parent structures, and writes canonical SMILES, InChI, InChIKey, molecular formula and exact molecular weight. Structure-quality statistics (parse success / failure counts) are emitted by this same script — there is no separate QC script. |
| `code_BioSurfDB/Model/model1/` | Task 1 — binary classification of a compound's **general biological function annotation status** (`biological_activity`, YES/NO) using Morgan fingerprints + Random Forest / XGBoost, with 5-fold stratified CV and bit-importance analysis. |
| `code_BioSurfDB/Model/model2/` | Task 2 — binary classification of a compound's **antimicrobial function annotation status** (`antimicrobial_activity`, YES/NO) using the same pipeline. |

Each module is self-contained and ships its own `environment.yml`,
`environment_full.yml` (a full locked snapshot for exact reproduction),
`requirements.txt` and `README_en.md` / `README_zh.md` with full setup and usage instructions.

## Quick start

The recommended setup is per-module — the environments are pre-configured with fast
mirrors and pin all versions. For example, to set up the RDKit module:

```bash
cd code_BioSurfDB/RDkit
conda env create -f environment.yml
conda activate surfdb
python scripts/deal_RDkit.py
```

> **Note on data:** each module's required input Excel files are included in this
> repository under that module's `input/` folder, so the scripts can run directly.
> Generated results are written to each module's `output/` folder (also tracked).
> Only Python build artifacts, Office lock files and logs are excluded — see
> `.gitignore`.

## Environments

Two environment files ship per module:
- `environment.yml` — a minimal, fast-to-build definition with all versions pinned and
  Tsinghua mirrors pre-configured (recommended for most users).
- `environment_full.yml` — a `conda env export` snapshot that pins **every** transitive
  dependency and build string. Use it only when you need bit-identical reproducibility
  of the environment used at publication time. Note it is platform-specific (captured on
  Windows).

## License

Released under the MIT License — see [LICENSE](LICENSE).
