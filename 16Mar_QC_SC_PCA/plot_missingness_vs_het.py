import pandas as pd
import matplotlib.pyplot as plt

# Define file paths
het_file = "/rds/general/user/wt1224/home/IG-heavy-chain-/16Mar_QC_SC_PCA/qc_step4_het.het"
miss_file = "/rds/general/user/wt1224/home/IG-heavy-chain-/16Mar_QC_SC_PCA/qc_step1_missing.smiss"

# Load heterozygosity data
het_df = pd.read_csv(het_file, sep='\s+')

# Load missingness data
miss_df = pd.read_csv(miss_file, sep='\s+')

# Merge datasets on FID and IID
merged_df = pd.merge(het_df, miss_df, on=["#FID", "IID"], how="inner")

# Print dataset shape to confirm merge success
print(f"Merged dataset shape: {merged_df.shape}")

# Create scatter plot
plt.figure(figsize=(8,6))
plt.scatter(merged_df["F_MISS"], merged_df["F"], alpha=0.5, edgecolors="k")
plt.xlabel("Missingness (F_MISS)")
plt.ylabel("Heterozygosity (F)")
plt.title("Missingness vs. Heterozygosity")
plt.grid(True)

# Save the plot
plot_path = "/rds/general/user/wt1224/home/IG-heavy-chain-/16Mar_QC_SC_PCA/missingness_vs_heterozygosity.png"
plt.savefig(plot_path, dpi=300)
print(f"Plot saved to: {plot_path}")