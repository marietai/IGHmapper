


import pandas as pd
import matplotlib.pyplot as plt

# Load missingness data
missing_data = pd.read_csv("/rds/general/user/wt1224/home/IG-heavy-chain-/my_genotype_data/final_missingness.smiss", sep="\s+")

# Calculate heterozygosity (1 - F_MISS)
missing_data["Heterozygosity"] = 1 - missing_data["F_MISS"]

# Plot
plt.figure(figsize=(8,6))
plt.scatter(missing_data["F_MISS"], missing_data["Heterozygosity"], alpha=0.5, color="blue")
plt.xlabel("Missingness (F_MISS)")
plt.ylabel("Heterozygosity")
plt.title("Missingness vs Heterozygosity")
plt.grid(True)

# Save plot
plt.savefig("missing_vs_het_plot.png")
print("Plot saved as missing_vs_het_plot.png")

