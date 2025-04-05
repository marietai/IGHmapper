import pandas as pd

# Load the map file and standardize chromosome format
map_file = r"C:\Users\marietai\Desktop\GoalPost\SNP\37_to_igh\igh_to_hg19_formatted.txt"
map_df = pd.read_csv(map_file, sep="\t", header=None, names=["Chromosome", "GRCh37_Pos", "GRCh37_Allele", "IGH_Pos", "IGH_Allele"])
map_df["Chromosome"] = map_df["Chromosome"].astype(str)

# Load the original .bim file
bim_file = r"C:\Users\marietai\Desktop\GoalPost\SNP\37_to_igh\GoalPost37_IGH_loci.bim"
bim_df = pd.read_csv(bim_file, sep="\t", header=None, names=["Chromosome", "SNP", "cM", "GRCh37_Pos", "Allele1", "Allele2"], dtype={"Chromosome": str})

# Ensure the chromosome in .bim file is set to '14' (You might want to verify this step)
bim_df['Chromosome'] = bim_df['Chromosome'].replace(['1'], '14')

# Merge on Chromosome and Position
merged_df = pd.merge(bim_df, map_df, how="left", left_on=["Chromosome", "GRCh37_Pos"], right_on=["Chromosome", "GRCh37_Pos"])

# Check for successful merging
if merged_df.empty:
    print("Error: No matching positions found between .bim file and map file.")
else:
    print("Merge successful, proceeding with update.")
    # Continue only if the merge was successful

    # Update positions and alleles where mapping is successful
    merged_df["GRCh37_Pos"] = merged_df["IGH_Pos"].fillna(merged_df["GRCh37_Pos"]).astype(int)
    merged_df["Allele1"] = merged_df["IGH_Allele"].fillna(merged_df["Allele1"])
    merged_df["Allele2"] = merged_df["IGH_Allele"].fillna(merged_df["Allele2"])
    
    # Drop unnecessary columns and save the updated .bim file
    output_bim_file = r"C:\Users\marietai\Desktop\GoalPost\SNP\37_to_igh\GoalPost37_IGH_reference.bim"
    merged_df[["Chromosome", "SNP", "cM", "GRCh37_Pos", "Allele1", "Allele2"]].to_csv(output_bim_file, sep="\t", index=False, header=False)
    
    print(f"Corrected .bim file saved to: {output_bim_file}")
