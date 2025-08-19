import pandas as pd
from scipy.stats import rankdata, norm

# === Step 1: Load the Excel file ===
input_path = "/mnt/c/Users/marietai/Desktop/GoalPost/phenotypic_analysis/Multiplex Serology GoalPost Data.xlsx"
df = pd.read_excel(input_path)

# === Step 2: Define rank normalization function ===
def rank_normalize(series):
    ranks = rankdata(series, method='average')
    n = len(ranks)
    return norm.ppf((ranks - 0.5) / n)

# === Step 3: Apply to antibody columns G to P ===
# These are columns 6 to 15 in 0-based indexing
antibody_cols = df.columns[6:16]
for col in antibody_cols:
    df[col + '_ranknorm'] = rank_normalize(df[col])

# === Step 4: Save output ===
output_tsv = "/mnt/c/Users/marietai/Desktop/GoalPost/phenotypic_analysis/Multiplex_Normalized.tsv"
df.to_csv(output_tsv, sep="\t", index=False)

print("✅ Rank normalization complete. Output saved.")
