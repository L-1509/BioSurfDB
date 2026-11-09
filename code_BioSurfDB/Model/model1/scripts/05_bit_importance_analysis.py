#!/usr/bin/env python
"""
Morgan bit importance analysis for RF and XGBoost.
Uses complete fingerprint files and 5-fold split files.
Writes permutation-importance bar plots and a per-bit summary table to OUT_DIR.
All file paths are absolute paths derived from the script location.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.metrics import roc_auc_score
import xgboost as xgb
import warnings
warnings.filterwarnings("ignore")

# ===================== Absolute path base =====================
BASE_DIR = Path(__file__).resolve().parent.parent

# ===================== Configuration =====================
N_SPLITS = 5
NBITS = 1024
RANDOM_STATE = 42
N_REPEATS = 5                # permutation importance repeats
TOP_N_PER_FOLD = 30          # top N bits per fold
MIN_STABLE_FOLDS = 3         # must appear in at least N folds

FP_DIR = BASE_DIR / "input" / "fingerprints_biological"
FOLDS_DIR = BASE_DIR / "input" / "cv_folds_biological"
TARGET_COL = "biological_activity"
OUT_DIR = BASE_DIR / "output" / "bit_importance_biological"

OUT_DIR.mkdir(parents=True, exist_ok=True)


# ===================== Load full fingerprints =====================
def load_all_fingerprints():
    pos_df = pd.read_csv(FP_DIR / "fingerprints_pos.csv")
    comp_df = pd.read_csv(FP_DIR / "fingerprints_comp.csv")

    df = pd.concat([pos_df, comp_df], ignore_index=True)
    df = df.set_index("sample_id")

    fp_cols = [c for c in df.columns if c.startswith("FP_")]
    if len(fp_cols) != NBITS:
        print(f"Warning: found {len(fp_cols)} FP columns, expected {NBITS}")

    return df, fp_cols


# ===================== Load fold data =====================
def load_fold(fold, fp_df, fp_cols):
    train_ids = pd.read_csv(FOLDS_DIR / f"fold_{fold}_train.csv")
    val_ids = pd.read_csv(FOLDS_DIR / f"fold_{fold}_val.csv")

    train_df = fp_df.loc[train_ids["sample_id"]]
    val_df = fp_df.loc[val_ids["sample_id"]]

    X_train = train_df[fp_cols].values.astype(np.uint8)
    y_train = train_df[TARGET_COL].values.astype(int)

    X_val = val_df[fp_cols].values.astype(np.uint8)
    y_val = val_df[TARGET_COL].values.astype(int)

    return X_train, y_train, X_val, y_val


# ===================== Models =====================
def make_rf():
    return RandomForestClassifier(
        n_estimators=100,
        random_state=RANDOM_STATE,
        n_jobs=-1,
        class_weight="balanced",
    )


def make_xgb():
    return xgb.XGBClassifier(
        n_estimators=100,
        random_state=RANDOM_STATE,
        eval_metric="logloss",
    )


# ===================== Plotting =====================
def plot_stable_bits(perm_mean, stability, model_name, save_path):
    stable_idx = np.where(stability >= MIN_STABLE_FOLDS)[0]
    if len(stable_idx) == 0:
        print(f"{model_name}: no stable bits found, try lowering MIN_STABLE_FOLDS")
        return

    stable_idx = stable_idx[np.argsort(perm_mean[stable_idx])[::-1]]
    top_idx = stable_idx[:20]

    values = perm_mean[top_idx]
    labels = [f"bit {i}" for i in top_idx]

    plt.figure(figsize=(8, 6))
    plt.barh(range(len(top_idx)), values[::-1], color="#1f77b4")
    plt.yticks(range(len(top_idx)), labels[::-1])
    plt.xlabel("Mean permutation importance (ROC-AUC decrease)")
    plt.title(f"Stable Morgan bits across folds: {model_name}")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches="tight", transparent=True)
    plt.close()
    print(f"Saved {save_path}")


# ===================== Main =====================
def main():
    print("=" * 70)
    print("Bit importance analysis - Task: biological")
    print(f"Base directory: {BASE_DIR}")
    print("=" * 70)

    fp_df, fp_cols = load_all_fingerprints()
    print(f"Total samples: {len(fp_df)}, fingerprint bits: {len(fp_cols)}")

    rf_perm_all = []
    xgb_perm_all = []
    rf_top_sets = []
    xgb_top_sets = []

    for fold in range(1, N_SPLITS + 1):
        print(f"\n===== Fold {fold} =====")
        X_train, y_train, X_val, y_val = load_fold(fold, fp_df, fp_cols)
        print(
            f"Train: {len(X_train)}, Val: {len(X_val)} "
            f"(Pos ratio: {y_val.sum() / len(y_val):.3f})"
        )

        # ---------- Random Forest ----------
        rf = make_rf()
        rf.fit(X_train, y_train)
        rf_auc = roc_auc_score(y_val, rf.predict_proba(X_val)[:, 1])
        rf_perm = permutation_importance(
            rf, X_val, y_val,
            scoring="roc_auc",
            n_repeats=N_REPEATS,
            random_state=RANDOM_STATE,
            n_jobs=-1,
        )
        rf_perm_mean = rf_perm.importances_mean
        rf_perm_all.append(rf_perm_mean)
        rf_top = np.argsort(rf_perm_mean)[::-1][:TOP_N_PER_FOLD]
        rf_top_sets.append(set(rf_top))
        print(f"RF  val AUC: {rf_auc:.4f}")

        # ---------- XGBoost ----------
        xgb_model = make_xgb()
        xgb_model.fit(X_train, y_train)
        xgb_auc = roc_auc_score(y_val, xgb_model.predict_proba(X_val)[:, 1])
        xgb_perm = permutation_importance(
            xgb_model, X_val, y_val,
            scoring="roc_auc",
            n_repeats=N_REPEATS,
            random_state=RANDOM_STATE,
            n_jobs=-1,
        )
        xgb_perm_mean = xgb_perm.importances_mean
        xgb_perm_all.append(xgb_perm_mean)
        xgb_top = np.argsort(xgb_perm_mean)[::-1][:TOP_N_PER_FOLD]
        xgb_top_sets.append(set(xgb_top))
        print(f"XGB val AUC: {xgb_auc:.4f}")

    # ===================== Aggregate across folds =====================
    rf_perm_all = np.array(rf_perm_all)     # (5, 1024)
    xgb_perm_all = np.array(xgb_perm_all)

    rf_perm_mean = rf_perm_all.mean(axis=0)
    xgb_perm_mean = xgb_perm_all.mean(axis=0)

    rf_stability = np.zeros(NBITS, dtype=int)
    for s in rf_top_sets:
        for b in s:
            rf_stability[b] += 1

    xgb_stability = np.zeros(NBITS, dtype=int)
    for s in xgb_top_sets:
        for b in s:
            xgb_stability[b] += 1

    # ===================== Plot =====================
    plot_stable_bits(
        rf_perm_mean, rf_stability, "Random Forest",
        OUT_DIR / "rf_stable_bits.png",
    )
    plot_stable_bits(
        xgb_perm_mean, xgb_stability, "XGBoost",
        OUT_DIR / "xgb_stable_bits.png",
    )

    # ===================== Save summary table =====================
    df_out = pd.DataFrame({
        "bit": np.arange(NBITS),
        "rf_perm_mean": rf_perm_mean,
        "rf_perm_std": rf_perm_all.std(axis=0),
        "rf_stability": rf_stability,
        "xgb_perm_mean": xgb_perm_mean,
        "xgb_perm_std": xgb_perm_all.std(axis=0),
        "xgb_stability": xgb_stability,
    })
    df_out.to_csv(
        OUT_DIR / "bit_importance_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )

    print(f"\nAll results saved to {OUT_DIR}")


if __name__ == "__main__":
    main()
