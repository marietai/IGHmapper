import pandas as pd

# Load the reference file (igh_to_hg38.txt)
map_file = r"C:\Users\marietai\Desktop\GoalPost\SNP\38_to_igh\igh_to_hg38.txt"
map_df = pd.read_csv(map_file, sep="\t", header=None, names=["Chromosome", "GRCh38_Pos", "GRCh38_Allele", "IGH_Pos", "IGH_Allele"])
map_df["Chromosome"] = map_df["Chromosome"].astype(str).str.replace("chr", "")

# Load the original .bim file for GRCh38 IGH locus
bim_file = r"C:\Users\marietai\Desktop\GoalPost\SNP\38_to_igh\GoalPost38_IGH_loci.bim"
bim_df = pd.read_csv(bim_file, sep="\t", header=None, names=["Chromosome", "SNP", "cM", "GRCh38_Pos", "Allele1", "Allele2"], dtype={"Chromosome": str})

# Ensure the chromosome in .bim file is '14'
bim_df['Chromosome'] = bim_df['Chromosome'].replace(['1'], '14')

# Merge on Chromosome and Position
merged_df = pd.merge(bim_df, map_df, how="left", left_on=["Chromosome", "GRCh38_Pos"], right_on=["Chromosome", "GRCh38_Pos"])

# Check for successful merging
if merged_df.empty:
    print("Error: No matching positions found between .bim file and map file.")
else:
    print("Merge successful, proceeding with update.")
    
    # Update positions and alleles where mapping is successful
    merged_df["GRCh38_Pos"] = merged_df["IGH_Pos"].fillna(merged_df["GRCh38_Pos"]).astype(int)
    merged_df["Allele1"] = merged_df["IGH_Allele"].fillna(merged_df["Allele1"])
    merged_df["Allele2"] = merged_df["IGH_Allele"].fillna(merged_df["Allele2"])
    
    # Drop unnecessary columns and save the updated .bim file
    output_bim_file = r"C:\Users\marietai\Desktop\GoalPost\SNP\38_to_igh\GoalPost38_IGH_reference.bim"
    merged_df[["Chromosome", "SNP", "cM", "GRCh38_Pos", "Allele1", "Allele2"]].to_csv(output_bim_file, sep="\t", index=False, header=False)
    
    print(f"Corrected .bim file saved to: {output_bim_file}")
