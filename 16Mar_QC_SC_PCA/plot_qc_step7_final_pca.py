import matplotlib.pyplot as plt
import numpy as np

# Load PCA results from the final QC step
pca_file = "/rds/general/user/wt1224/home/IG-heavy-chain-/16Mar_QC_SC_PCA/qc_step7_final_pca.eigenvec"
data = np.loadtxt(pca_file, dtype=str)

# Extract IDs and PCA components
ids = data[:, :2]  # First two columns contain IDs
pca_values = data[:, 2:].astype(float)  # Convert PCA values to float

# Plot PC1 vs PC2
plt.figure(figsize=(8,6))
plt.scatter(pca_values[:, 0], pca_values[:, 1], alpha=0.5, edgecolors="k")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("Final PCA: Population Structure")
plt.grid(True)

# Save the plot
output_plot = "/rds/general/user/wt1224/home/IG-heavy-chain-/16Mar_QC_SC_PCA/qc_step7_final_pca_plot.png"
plt.savefig(output_plot, dpi=300)
print(f"Final PCA plot saved at: {output_plot}")
