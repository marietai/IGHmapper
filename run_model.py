#!/usr/bin/env python3
"""
run_model.py — Score new samples with a trained IGH classifier and apply thresholds.

What it does
------------
- Loads a trained pipeline (.joblib) and its threshold JSON (from train_model.py).
- Reads a TSV of samples, aligns feature columns to the model, and computes
  P(label==0) = probability of being an "accurate" SNP ("p_good").
- Applies either:
    (A) a single unique threshold (from the JSON), or
    (B) stratified thresholds: HIGH threshold for a high-risk stratum and ELSE threshold otherwise.
- Writes a TSV with probabilities and predicted labels.

Stratified rule (default)
-------------------------
HIGH stratum: (sv_del == 0) OR (sv_complex == 0) OR (igh_family in {1,3,4})
Thresholds:
  --thr-high (default 0.99)   applied to HIGH stratum
  --thr-else (default: use JSON chosen threshold) applied to ELSE stratum

Examples
--------
# Unique threshold (from JSON)
python run_model.py \
  --model model_11/RandomForest_model.joblib \
  --threshold-json model_11/RandomForest_threshold.json \
  --input testing_sets/testing_features_allsamples.tsv \
  --output model_11/model_11_prediction_unique_thr.tsv

# Stratified thresholds
python run_model.py \
  --model model_11/RandomForest_model.joblib \
  --threshold-json model_11/RandomForest_threshold.json \
  --input testing_sets/testing_features_allsamples.tsv \
  --output model_11/model_11_prediction_stratified_thr.tsv \
  --use-stratified \
  --thr-high 0.99
"""

from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
from typing import Optional

import joblib
import numpy as np
import pandas as pd


# ----------------------------
# CLI
# ----------------------------
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Score samples and apply unique/stratified thresholds.")
    p.add_argument("--model", required=True, help="Path to trained model pipeline (.joblib).")
    p.add_argument("--threshold-json", required=True, help="Threshold JSON saved during training.")
    p.add_argument("--input", required=True, help="Input TSV with features.")
    p.add_argument("--output", required=True, help="Output TSV with probabilities and predictions.")
    p.add_argument(
        "--drop-cols",
        nargs="*",
        default=[
            "igh_sample","ref","alt","gt","position",
            "abs_baf_diff","abs_gc_diff","abs_log2r","log_dist_to_sv",
            "label","lr_REF","lr_ALT","lr_GT"
        ],
        help="Columns to drop from the input before scoring (should mirror training).",
    )
    p.add_argument("--use-stratified", action="store_true", help="Use stratified thresholds instead of unique.")
    p.add_argument("--thr-high", type=float, default=0.99, help="Threshold for HIGH stratum.")
    p.add_argument(
        "--thr-else",
        type=float,
        default=None,
        help="Threshold for ELSE stratum. Default: use the JSON chosen threshold.",
    )
    # Advanced: allow overriding the high-stratum definition via simple flags
    p.add_argument("--high-on-del", action="store_true", default=True, help="HIGH if sv_del == 0 (default: on).")
    p.add_argument("--high-on-complex", action="store_true", default=True, help="HIGH if sv_complex == 0 (default: on).")
    p.add_argument("--high-on-families", nargs="*", type=int, default=[1,3,4],
                   help="IGH families marking HIGH stratum (default: 1 3 4).")
    return p.parse_args()


# ----------------------------
# Utilities
# ----------------------------
def _align_columns_like_model(X: pd.DataFrame, model) -> pd.DataFrame:
    """
    Reorder and/or add missing columns (filled with NaN -> imputer handles) to match the model input.
    """
    feature_names: Optional[list[str]] = None

    # Try pipeline final step for feature_names_in_
    try:
        last = model.named_steps[list(model.named_steps.keys())[-1]]
        if hasattr(last, "feature_names_in_"):
            feature_names = list(last.feature_names_in_)
    except Exception:
        pass

    # Try top-level (some sklearn versions store at pipeline level)
    if feature_names is None and hasattr(model, "feature_names_in_"):
        feature_names = list(model.feature_names_in_)

    if feature_names is not None:
        missing = [c for c in feature_names if c not in X.columns]
        for c in missing:
            X[c] = np.nan
        X = X[feature_names]
    else:
        # Fallback: stable sorted order (must match training choice)
        X = X.reindex(sorted(X.columns), axis=1)

    return X


def _get_pos_class_index(model, class_index_from_json: Optional[int]) -> int:
    """
    Return the column index of the positive class (label==0) for predict_proba.
    Prefer stored JSON index; otherwise infer from estimator classes_.
    """
    if class_index_from_json is not None:
        return int(class_index_from_json)

    # Try pipeline final estimator
    if hasattr(model, "named_steps"):
        last = model.named_steps[list(model.named_steps.keys())[-1]]
        if hasattr(last, "classes_"):
            return int(np.where(last.classes_ == 0)[0][0])

    # Try top-level
    if hasattr(model, "classes_"):
        return int(np.where(model.classes_ == 0)[0][0])

    raise RuntimeError("Could not locate classes_ on the model; please supply 'class_index_used' in the JSON.")


# ----------------------------
# Main
# ----------------------------
def main():
    args = parse_args()

    # Load model and threshold record
    model = joblib.load(args.model)
    thr_info = json.loads(Path(args.threshold_json).read_text())
    thr_unique = float(thr_info["chosen_threshold"])
    pos_col_idx = _get_pos_class_index(model, thr_info.get("class_index_used"))

    # Load input (do NOT drop rows; let imputer handle missing)
    df = pd.read_csv(args.input, sep="\t")
    drop_cols = [c for c in args.drop_cols if c in df.columns]
    X = df.drop(columns=drop_cols, errors="ignore")

    # Align columns to model's expected inputs
    X = _align_columns_like_model(X, model)

    # Probabilities of positive class (label==0 = "accurate")
    proba_pos = model.predict_proba(X)[:, pos_col_idx]

    # Output frame with probability
    out = df.copy()
    out["p_good"] = proba_pos  # consistent naming for downstream logic

    # Threshold selection mode
    if not args.use_stratified:
        # Unique threshold mode
        preds = (out["p_good"] >= thr_unique).astype(int)  # 1 = predicted accurate (label 0)
        mode_desc = f"UNIQUE threshold: {thr_unique:.3f} (target PPV={thr_info.get('target_ppv')})"
        n_high = n_else = None
    else:
        # Stratified thresholds
        req = ["sv_del", "sv_complex", "igh_family"]
        missing = [c for c in req if c not in out.columns]
        if missing:
            raise ValueError(f"Missing required columns for stratified thresholding: {missing}")

        sv_del = pd.to_numeric(out["sv_del"], errors="coerce")
        sv_complex = pd.to_numeric(out["sv_complex"], errors="coerce")
        igh_family = pd.to_numeric(out["igh_family"], errors="coerce")

        high_mask = np.zeros(len(out), dtype=bool)
        if args.high_on_del:
            high_mask |= (sv_del == 0)
        if args.high_on_complex:
            high_mask |= (sv_complex == 0)
        if args.high_on_families:
            high_mask |= igh_family.isin(args.high_on_families).to_numpy()

        thr_else = float(args.thr_else) if args.thr_else is not None else thr_unique
        preds = np.where(high_mask, (out["p_good"] >= args.thr_high), (out["p_good"] >= thr_else)).astype(int)

        n_high = int(high_mask.sum())
        n_else = int(len(out) - n_high)
        mode_desc = (
            f"STRATIFIED thresholds: HIGH={args.thr_high:.3f} when "
            f"(sv_del==0 OR sv_complex==0 OR igh_family in {set(args.high_on_families)}) [n={n_high}], "
            f"ELSE={thr_else:.3f} [n={n_else}] "
            f"(base unique thr={thr_unique:.3f}, target PPV={thr_info.get('target_ppv')})"
        )

    print(f"Using: {mode_desc}")

    # Save results
    out["probability"] = out["p_good"]
    out["predicted_label"] = preds  # 1 = predicted accurate (label 0), 0 = predicted inaccurate (label 1)
    out_dir = Path(args.output).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False)
    print(f"✅ Predictions saved to: {args.output}")

    # Summary
    n_pos = int(preds.sum())
    rate = n_pos / len(out) if len(out) else 0.0
    print(f"Predicted 'accurate' (label 0): {n_pos} / {len(out)} ({rate:.2%})")
    if args.use_stratified and n_high is not None:
        kept_high = int(((out["p_good"] >= args.thr_high) & (preds == 1)).sum())
        kept_else = int(((out["p_good"] >= (args.thr_else if args.thr_else is not None else thr_unique)) & (preds == 1) & ~((out["sv_del"]==0) | (out["sv_complex"]==0) | (out["igh_family"].isin(args.high_on_families)))).sum())
        print(f"  Kept in HIGH stratum: {kept_high} / {n_high}")
        print(f"  Kept in ELSE stratum: {kept_else} / {n_else}")


if __name__ == "__main__":
    main()

