#!/usr/bin/env python
"""
Stratified 5-fold cross-validation split based on the fingerprint files.
Only metadata and sample_id are saved for each fold.
All file paths are absolute paths derived from the script location.
"""

import os
from pathlib import Path

import pandas as pd
from sklearn.model_selection import StratifiedKFold

# ===================== Absolute path base =====================
BASE_DIR = Path(__file__).resolve().parent.parent

# ===================== Configuration =====================
FP_DIR = BASE_DIR / "input" / "fingerprints_biological"
OUTPUT_DIR = BASE_DIR / "input" / "cv_folds_biological"

RANDOM_STATE = 42
N_SPLITS = 5
STRATIFY_COL = "biological_activity"

META_COLS = [
    "sample_id",
    "Database_ID",
    "SMILES_Parent_RDKit",
    "biological_activity",
    "antimicrobial_activity",
    "source",
]

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ===================== Load fingerprint files =====================
print("Reading fingerprint files...")
print(f"Base directory: {BASE_DIR}")

pos_df = pd.read_csv(FP_DIR / "fingerprints_pos.csv")
comp_df = pd.read_csv(FP_DIR / "fingerprints_comp.csv")

df_all = pd.concat([pos_df, comp_df], ignore_index=True)
df_all[STRATIFY_COL] = df_all[STRATIFY_COL].astype(int)

print(f"Total samples: {len(df_all)}")
print(
    f"Pos: {df_all[STRATIFY_COL].sum()}, "
    f"Comp: {len(df_all) - df_all[STRATIFY_COL].sum()}"
)

# ===================== Stratified 5-fold split =====================
y = df_all[STRATIFY_COL].values

skf = StratifiedKFold(
    n_splits=N_SPLITS,
    shuffle=True,
    random_state=RANDOM_STATE,
)

for fold_num, (train_index, val_index) in enumerate(
    skf.split(df_all, y), start=1
):
    train_df = df_all.iloc[train_index][META_COLS].copy()
    val_df = df_all.iloc[val_index][META_COLS].copy()

    train_path = OUTPUT_DIR / f"fold_{fold_num}_train.csv"
    val_path = OUTPUT_DIR / f"fold_{fold_num}_val.csv"

    train_df.to_csv(train_path, index=False)
    val_df.to_csv(val_path, index=False)

    train_pos = train_df[STRATIFY_COL].sum()
    train_comp = len(train_df) - train_pos
    val_pos = val_df[STRATIFY_COL].sum()
    val_comp = len(val_df) - val_pos

    print(
        f"Fold {fold_num}: "
        f"Train {len(train_df)} (pos:{train_pos}, comp:{train_comp}), "
        f"Val {len(val_df)} (pos:{val_pos}, comp:{val_comp})"
    )

print(f"\nSplit complete! Files saved to: {OUTPUT_DIR}/")
