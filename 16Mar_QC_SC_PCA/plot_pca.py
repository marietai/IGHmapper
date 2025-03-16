import matplotlib.pyplot as plt
import numpy as np

# Load PCA results
pca_file = "/rds/general/user/wt1224/home/IG-heavy-chain-/16Mar_QC_SC_PCA/qc_step6_pca.eigenvec"
data = np.loadtxt(pca_file, dtype=str)

# Extract IDs and PCA components
ids = data[:, :2]
pca_values = data[:, 2:].astype(float)

# Plot PC1 vs PC2
plt.figure(figsize=(8,6))
plt.scatter(pca_values[:, 0], pca_values[:, 1], alpha=0.5)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("PCA: Population Structure")
plt.savefig("/rds/general/user/wt1224/home/IG-heavy-chain-/16Mar_QC_SC_PCA/qc_step6_pca_plot.png")

