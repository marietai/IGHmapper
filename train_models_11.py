import pandas as pd
import numpy as np
import gc
import mlflow
import mlflow.sklearn
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import joblib
warnings.filterwarnings("ignore")
import json
from datetime import datetime

# === Feature name mapping for plotting ===
FEATURE_NAME_MAP = {
    "IGH_Pos": "IGH custom reference coordinate",
    "gc_content": "GC content (%)",
    "baf": "BAF",
    "log2r": "Log2 ratio",
    "sv_del": "SNP overlaps deletion SV",
    "sv_dup": "SNP overlaps duplication SV",
    "sv_msv": "SNP overlaps multi-allelic SV",
    "sv_complex": "SNP overlaps complex SV",
    "sv_count": "SV overlap count",
    "in_igh_gene": "SNP within IGH gene",
    "igh_family": "IGH gene family",
    "igh_uniqueness_score": "IGH gene positional redundancy score",
    "near_igh_gene": "SNP near IGH gene (±100 bp)",
    "confidence": "Confidence score",
    "true_dist_to_sv": "Distance to nearest SV (kb)"
}

from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.metrics import (
    precision_score, recall_score,
    average_precision_score, roc_auc_score,
    roc_curve, precision_recall_curve
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from matplotlib.colors import ListedColormap
from xgboost import XGBClassifier
import itertools
import os

plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.color": "black",
    "ytick.color": "black",
    "axes.edgecolor": "black",
    "axes.labelcolor": "black"
})

# === Load training and test sets explicitly ===
train_path = "/mnt/c/Users/marietai/Desktop/GoalPost/ML/data/training_features.127labeled.tsv"
test_path = "/mnt/c/Users/marietai/Desktop/GoalPost/ML/data/testing_features.48labeled.tsv"

# Columns to drop after handling missing data
drop_cols = ["igh_sample", "ref", "alt", "gt", "position","log_dist_to_sv", "abs_log2r", "abs_baf_diff", "abs_gc_diff", "lr_REF", "lr_ALT", "lr_GT", "label"]

# === Load and preprocess training set ===
df_train = pd.read_csv(train_path, sep="\t")
df_train = df_train.dropna()  # Drop rows with any missing values
y_train = df_train["label"]
X_train = df_train.drop(columns=drop_cols)

# === Load and preprocess test set ===
df_test = pd.read_csv(test_path, sep="\t")
df_test = df_test.dropna()  # Drop rows with any missing values
y_test = df_test["label"]
X_test = df_test.drop(columns=drop_cols)


# === Encode categorical features consistently between train and test ===
for col in X_train.select_dtypes(include=["object", "category"]).columns:
    le = LabelEncoder()
    all_values = pd.concat([X_train[col], X_test[col]], axis=0).astype(str)
    le.fit(all_values)
    X_train[col] = le.transform(X_train[col].astype(str))
    X_test[col] = le.transform(X_test[col].astype(str))

# === Model dictionary ===
models = {
    "RandomForest": RandomForestClassifier(n_estimators=500, random_state=42),
    "XGBoost_1Tree": XGBClassifier(n_estimators=500, use_label_encoder=False, eval_metric='logloss', random_state=42)
}

# === Create output folder ===
os.makedirs("model_11", exist_ok=True)

# === Training loop ===
for name, estimator in models.items():
    mlflow.set_experiment("IGH_11")
    model = make_pipeline(SimpleImputer(strategy="median"), estimator)  # ✅ define before fitting

    try:
        if mlflow.active_run():
            mlflow.end_run()
    except Exception as e:
        print(f"Warning: Could not end previous MLflow run. Reason: {e}")

    with mlflow.start_run(run_name=f"{name}_Split"):

        model.fit(X_train, y_train) 
        
        # === Save model ===
        model_filename = f"model_11/{name}_model.joblib"
        joblib.dump(model, model_filename)
        mlflow.log_artifact(model_filename)
        print(f"✅ Saved model to: {model_filename}")

        # === Positive class is label 0 (accurate) ===
        # Get column index for class 0 safely (pipeline or estimator)
        cls = model.classes_ if hasattr(model, "classes_") else model.named_steps[list(model.named_steps.keys())[-1]].classes_
        pos_idx = int(np.where(cls == 0)[0][0])

        # === Probabilities for POSITIVE class (label==0) on TRAIN & TEST ===
        proba_pos_tr = model.predict_proba(X_train)[:, pos_idx]
        proba_pos_te = model.predict_proba(X_test)[:,  pos_idx]

        # y_true remapped so that label==0 becomes 1 (positive) for metrics that expect positive=1
        y_pos_tr = (y_train == 0).astype(int)
        y_pos_te = (y_test  == 0).astype(int)

        # === AUCs (train & test) ===
        ap_tr   = average_precision_score(y_pos_tr, proba_pos_tr)   # PR AUC (train)
        rocau_tr = roc_auc_score(y_pos_tr, proba_pos_tr)            # ROC AUC (train)
        ap_te   = average_precision_score(y_pos_te, proba_pos_te)   # PR AUC (test)
        rocau_te = roc_auc_score(y_pos_te, proba_pos_te)            # ROC AUC (test)

        # === Choose a precision-first test threshold (e.g., PPV >= 0.99) ===
        prec_curve_te, rec_curve_te, thr_curve_te = precision_recall_curve(y_pos_te, proba_pos_te)  # prec/rec: n+1, thr: n
        prec_te = prec_curve_te[1:]
        rec_te  = rec_curve_te[1:]
        thr_te  = thr_curve_te

        target_ppv = 0.99
        mask = prec_te >= target_ppv
        if mask.any():
            idx = np.argmax(rec_te[mask])               # among PPV>=target, maximize recall
            chosen_thr = float(thr_te[mask][idx])
            op_idx = int(np.where(mask)[0][idx])        # index in prec/rec (aligned to thr)
        else:
            beta = 0.5
            fbeta = (1 + beta**2) * (prec_te * rec_te) / (beta**2 * prec_te + rec_te + 1e-12)
            idx = int(np.nanargmax(fbeta))
            chosen_thr = float(thr_te[idx])
            op_idx = idx

        # Predictions at the chosen test threshold (1 = predicted accurate/label 0)
        preds_te = (proba_pos_te >= chosen_thr).astype(int)

        # Precision/Recall at threshold (test)
        ppv_te = precision_score(y_pos_te, preds_te, zero_division=0)
        tpr_te = recall_score(y_pos_te, preds_te, zero_division=0)

        # ... after chosen_thr is computed and before plotting/logging:
        thr_record = {
            "model_name": name,
            "target_ppv": float(target_ppv),
            "chosen_threshold": float(chosen_thr),
            "pos_label": 0,                       # positive = label 0 (accurate)
            "class_index_used": int(pos_idx),     # which column in predict_proba
            "timestamp": datetime.now().isoformat()
        }

        thr_path = f"model_11/{name}_threshold.json"
        with open(thr_path, "w") as f:
            json.dump(thr_record, f, indent=2)
        print(f"💾 Saved threshold to: {thr_path}")
        mlflow.log_artifact(thr_path)  # optional: keep with the run


        # === Log AUCs (train & test) and thresholded test metrics to MLflow ===
        mlflow.log_metrics({
            "pr_auc_pos0_train":  ap_tr,
            "roc_auc_pos0_train": rocau_tr,
            "pr_auc_pos0_test":   ap_te,
            "roc_auc_pos0_test":  rocau_te,
            "precision_at_thr_pos0_test": ppv_te,
            "recall_at_thr_pos0_test":    tpr_te,
            "chosen_threshold_pos0_test": chosen_thr
        })

        # (Optional) confusion matrix counts (test)
        from sklearn.metrics import confusion_matrix
        tn, fp, fn, tp = confusion_matrix(y_pos_te, preds_te, labels=[0, 1]).ravel()
        mlflow.log_metrics({"TP_pos0_test": int(tp), "FP_pos0_test": int(fp), "TN_pos0_test": int(tn), "FN_pos0_test": int(fn)})

        # === ROC curves (train & test; pos=label0) ===
        fpr_tr, tpr_tr, _ = roc_curve(y_pos_tr, proba_pos_tr)
        fpr_te, tpr_te_curve, _ = roc_curve(y_pos_te, proba_pos_te)

        plt.figure()
        plt.plot(fpr_tr, tpr_tr,        label=f"Train AUC = {rocau_tr:.2f}", color="black")
        plt.plot(fpr_te, tpr_te_curve,  label=f"Test  AUC = {rocau_te:.2f}", color="gray")
        plt.plot([0,1],[0,1],"--", color="lightgray")
        plt.xlabel("False positive rate"); plt.ylabel("True positive rate")
        plt.legend(loc="lower right", frameon=False)
        plt.tight_layout(); plt.savefig(f"model_11/{name}_ROC_curve_train_test.png", dpi=300); plt.close()

        # === PR curves (train & test) ===
        prec_tr, rec_tr, _ = precision_recall_curve(y_pos_tr, proba_pos_tr)
        prec_te, rec_te, _ = precision_recall_curve(y_pos_te, proba_pos_te)

        plt.figure()
        plt.plot(rec_tr, prec_tr, label=f"Train AP = {ap_tr:.2f}", color="black")
        plt.plot(rec_te, prec_te, label=f"Test  AP = {ap_te:.2f}", color="gray")
        # Mark test operating point
        plt.scatter(rec_curve_te[op_idx + 1], prec_curve_te[op_idx + 1], s=35)
        plt.text(rec_curve_te[op_idx + 1], prec_curve_te[op_idx + 1],
                f"thr={chosen_thr:.2f}\nPPV={ppv_te:.2f}, Rec={tpr_te:.2f}",
                fontsize=8, ha="right", va="top")
        plt.xlabel("Recall"); plt.ylabel("Precision")
        plt.legend(loc="lower left", frameon=False)
        plt.tight_layout()
        plt.savefig(f"model_11/{name}_PR_curve_train_test.png", dpi=300)
        plt.close()
        
           # === Feature Importance Plot with refined names ===
        try:
            final_model = model.named_steps[list(model.named_steps.keys())[-1]]
            if hasattr(final_model, "feature_importances_"):
                importances = final_model.feature_importances_
                feature_names = np.array([FEATURE_NAME_MAP.get(col, col) for col in X_train.columns])
                sorted_idx = np.argsort(importances)[::-1][:20]

                plt.figure(figsize=(6, 5))
                sns.barplot(
                    x=importances[sorted_idx],
                    y=feature_names[sorted_idx],
                    palette=["black"] + ["gray"] * 19
                )
                plt.xlabel("Feature Importance")
                #plt.title("Top 18 Features")
                plt.tight_layout()
                plt.savefig(f"model_11/{name}_feature_importance.png", dpi=300)
                plt.close()
        except Exception as e:
            print(f"Could not plot feature importance for {name}: {e}")

