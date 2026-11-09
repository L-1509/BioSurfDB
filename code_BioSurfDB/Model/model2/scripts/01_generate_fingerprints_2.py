#!/usr/bin/env python
"""
Task 2: Antimicrobial Activity.
Generate complete Morgan fingerprints.
Output two separate files for positive and comparison samples.
SMILES source column: SMILES_Parent_RDKit
All file paths are absolute paths derived from the script location.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs

# ===================== Absolute path base =====================
BASE_DIR = Path(__file__).resolve().parent.parent

# ===================== Configuration =====================
POS_FILE = BASE_DIR / "input" / "antimicrobial_pos.xlsx"
COMP_FILE = BASE_DIR / "input" / "antimicrobial_comp.xlsx"
OUTPUT_DIR = BASE_DIR / "input" / "fingerprints_antimicrobial"

RADIUS = 2
NBITS = 1024

SMILES_COL = "SMILES_Parent_RDKit"

KEEP_COLS = [
    "Database_ID",
    "SMILES_Parent_RDKit",
    "antimicrobial_activity",
]

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(
    radius=RADIUS,
    fpSize=NBITS,
)


# ===================== Utility functions =====================
def map_label(val):
    """Map various label formats to 1 or 0."""
    if pd.isna(val):
        return 0
    if isinstance(val, bool):
        return int(val)
    if isinstance(val, (int, float)):
        return 1 if val != 0 else 0
    if isinstance(val, str):
        v = val.strip().lower()
        return 1 if v in ["yes", "true", "1"] else 0
    return 0


def fp_to_numpy(fp):
    """Convert RDKit ExplicitBitVect to numpy uint8 array."""
    arr = np.zeros((NBITS,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def process_one_file(path, source):
    print(f"Reading {source}: {path}")
    df = pd.read_excel(path, sheet_name=0)
    df = df[KEEP_COLS].copy()

    df["antimicrobial_activity"] = df["antimicrobial_activity"].apply(map_label)

    df.dropna(subset=[SMILES_COL], inplace=True)
    df = df.reset_index(drop=True)

    df["source"] = source
    df["sample_id"] = [f"{source}_{i:06d}" for i in range(len(df))]

    fp_rows = []
    valid_idx = []

    for i, smi in enumerate(df[SMILES_COL].astype(str)):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"  [skip] invalid SMILES at row {i}: {smi}")
            continue
        fp = MORGAN_GEN.GetFingerprint(mol)
        fp_rows.append(fp_to_numpy(fp))
        valid_idx.append(i)

    if len(fp_rows) == 0:
        raise RuntimeError(
            f"No valid SMILES found in '{source}' file: {path}. "
            f"Check the SMILES column '{SMILES_COL}' and the input file."
        )
    df = df.iloc[valid_idx].reset_index(drop=True)
    fps = np.vstack(fp_rows)

    fp_cols = [f"FP_{i:04d}" for i in range(NBITS)]
    fp_df = pd.DataFrame(fps, columns=fp_cols)

    out_df = pd.concat([df, fp_df], axis=1)

    out_path = OUTPUT_DIR / f"fingerprints_{source}.csv"
    out_df.to_csv(out_path, index=False)
    print(f"{source}: {len(out_df)} samples saved to {out_path}")
    return out_df


if __name__ == "__main__":
    print("=" * 70)
    print("Task 2: Antimicrobial Activity - Generate Morgan fingerprints")
    print(f"SMILES source column: {SMILES_COL}")
    print(f"Base directory: {BASE_DIR}")
    print("=" * 70)

    pos_df = process_one_file(POS_FILE, "pos")
    comp_df = process_one_file(COMP_FILE, "comp")

    print("\nDone.")
    print(f"Output dir: {OUTPUT_DIR}")
