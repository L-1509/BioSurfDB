#!/usr/bin/env python
"""
Task 2: Antimicrobial Activity.
Random Forest 5-fold cross-validation.
Fingerprints are read from fingerprints_pos.csv / fingerprints_comp.csv.
All file paths are absolute paths derived from the script location.
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    matthews_corrcoef,
    balanced_accuracy_score,
    roc_curve,
)

import warnings
warnings.filterwarnings("ignore")

# ===================== Absolute path base =====================
BASE_DIR = Path(__file__).resolve().parent.parent

# ===================== Configuration =====================
FP_DIR = BASE_DIR / "input" / "fingerprints_antimicrobial"
FOLDS_DIR = BASE_DIR / "input" / "cv_folds_antimicrobial"
OUTPUT_DIR = BASE_DIR / "output" / "antimicrobial"

N_SPLITS = 5
RANDOM_STATE = 42
TARGET_COL = "antimicrobial_activity"

(OUTPUT_DIR / "figures").mkdir(parents=True, exist_ok=True)


# ===================== Data loading =====================
def load_all_fingerprints():
    pos_df = pd.read_csv(FP_DIR / "fingerprints_pos.csv")
    comp_df = pd.read_csv(FP_DIR / "fingerprints_comp.csv")

    df = pd.concat([pos_df, comp_df], ignore_index=True)
    fp_cols = [c for c in df.columns if c.startswith("FP_")]

    df = df.set_index("sample_id")
    return df, fp_cols


def load_fold_data(fold_num, fp_df, fp_cols):
    train_ids = pd.read_csv(FOLDS_DIR / f"fold_{fold_num}_train.csv")
    val_ids = pd.read_csv(FOLDS_DIR / f"fold_{fold_num}_val.csv")

    train_df = fp_df.loc[train_ids["sample_id"]].reset_index()
    val_df = fp_df.loc[val_ids["sample_id"]].reset_index()

    X_train = train_df[fp_cols].values.astype(np.uint8)
    y_train = train_df[TARGET_COL].values.astype(int)

    X_val = val_df[fp_cols].values.astype(np.uint8)
    y_val = val_df[TARGET_COL].values.astype(int)

    return X_train, y_train, X_val, y_val


def evaluate_model(model, X_train, y_train, X_val, y_val):
    model.fit(X_train, y_train)

    y_pred = model.predict(X_val)
    y_prob = model.predict_proba(X_val)[:, 1]

    metrics = {
        "accuracy": accuracy_score(y_val, y_pred),
        "precision": precision_score(y_val, y_pred, zero_division=0),
        "recall": recall_score(y_val, y_pred, zero_division=0),
        "f1": f1_score(y_val, y_pred, zero_division=0),
        "roc_auc": roc_auc_score(y_val, y_prob),
        "mcc": matthews_corrcoef(y_val, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_val, y_pred),
    }
    return metrics, y_prob


def plot_roc_curves(all_fpr, all_tpr, all_auc, save_path):
    plt.figure(figsize=(8, 7))
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
    linestyles = ["-", "--", "-.", ":", "-"]

    for i, (fpr, tpr, auc_val) in enumerate(zip(all_fpr, all_tpr, all_auc)):
        plt.plot(
            fpr, tpr,
            color=colors[i % len(colors)],
            linestyle=linestyles[i % len(linestyles)],
            linewidth=2,
            label=f"Fold {i + 1} (AUC = {auc_val:.3f})",
        )

    plt.plot([0, 1], [0, 1], "k--", linewidth=1,
             label="Random Classifier (AUC = 0.5)")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate", fontsize=12)
    plt.ylabel("True Positive Rate", fontsize=12)
    plt.title("Random Forest ROC - Antimicrobial 5-Fold CV", fontsize=14)
    plt.legend(loc="lower right", fontsize=10)
    plt.grid(alpha=0.3)
    plt.savefig(save_path, dpi=300, bbox_inches="tight", transparent=True)
    plt.close()
    print(f"ROC curves saved to {save_path}")


# ===================== Main =====================
if __name__ == "__main__":
    print("=" * 70)
    print(f"Task 2: Antimicrobial Activity - Random Forest 5-Fold CV")
    print(f"Base directory: {BASE_DIR}")
    print("=" * 70)

    fp_df, fp_cols = load_all_fingerprints()
    print(f"Fingerprint bits: {len(fp_cols)}")

    model = RandomForestClassifier(
        n_estimators=100,
        random_state=RANDOM_STATE,
        n_jobs=-1,
        class_weight="balanced",
    )

    fold_metrics = []
    all_fpr, all_tpr, all_auc = [], [], []

    for fold in range(1, N_SPLITS + 1):
        print(f"\n--- Fold {fold} ---")
        X_train, y_train, X_val, y_val = load_fold_data(fold, fp_df, fp_cols)
        print(
            f"Train: {len(X_train)}, Val: {len(X_val)} "
            f"(Positive ratio: {y_val.sum() / len(y_val):.3f})"
        )

        metrics, y_prob = evaluate_model(
            model, X_train, y_train, X_val, y_val
        )
        fold_metrics.append(metrics)

        print(
            f"Acc:{metrics['accuracy']:.3f} "
            f"Prec:{metrics['precision']:.3f} "
            f"Rec:{metrics['recall']:.3f} "
            f"F1:{metrics['f1']:.3f} "
            f"AUC:{metrics['roc_auc']:.3f} "
            f"MCC:{metrics['mcc']:.3f} "
            f"BAcc:{metrics['balanced_accuracy']:.3f}"
        )

        fpr, tpr, _ = roc_curve(y_val, y_prob)
        all_fpr.append(fpr)
        all_tpr.append(tpr)
        all_auc.append(metrics["roc_auc"])

    plot_roc_curves(
        all_fpr, all_tpr, all_auc,
        save_path=OUTPUT_DIR / "figures" / "rf_roc_curves_2.png",
    )

    df_metrics = pd.DataFrame(fold_metrics)
    mean = df_metrics.mean()
    std = df_metrics.std()

    print("\n" + "=" * 70)
    print("5-Fold Summary (Mean +/- Std)")
    print("=" * 70)
    for metric in df_metrics.columns:
        print(f"{metric:15}: {mean[metric]:.3f} +/- {std[metric]:.3f}")

    df_metrics.to_csv(OUTPUT_DIR / "rf_cv_results_2.csv", index=False)
    mean.to_frame(name="mean").T.to_csv(OUTPUT_DIR / "rf_summary_2.csv")

    print("\nResults saved.")
