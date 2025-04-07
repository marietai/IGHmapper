import os
import pandas as pd

def load_mapping_file(build_type):
    base_path = os.path.dirname(__file__)
    if build_type == "GRCh37":
        map_file = os.path.join(base_path, "igh_to_hg19.txt")
    elif build_type == "GRCh38":
        map_file = os.path.join(base_path, "igh_to_hg38.txt")
    else:
        raise ValueError("Invalid build type. Choose 'GRCh37' or 'GRCh38'.")
    
    map_df = pd.read_csv(map_file, sep="\t", header=None, names=["Chromosome", "GRCh38_Pos", "GRCh38_Allele", "IGH_Pos", "IGH_Allele"])
    map_df["Chromosome"] = map_df["Chromosome"].str.replace("chr", "")
    return map_df

def convert_to_custom_reference(bim_file, build_type, output_bim_file):
    bim_df = pd.read_csv(bim_file, sep="\t", header=None, names=["Chromosome", "SNP", "cM", "Position", "Allele1", "Allele2"], dtype={"Chromosome": str})

    # Convert Chromosome column to string type
    bim_df["Chromosome"] = bim_df["Chromosome"].astype(str)
    
    # Load the mapping file
    map_df = load_mapping_file(build_type)

    # Ensure Chromosome column is compatible
    map_df["Chromosome"] = map_df["Chromosome"].astype(str)

    # Merge dataframes
    merged_df = pd.merge(bim_df, map_df, how="left", left_on=["Chromosome", "Position"], right_on=["Chromosome", "GRCh38_Pos"])

    # Replace positions and alleles in .bim file with IGH_Pos and IGH_Allele where mapping exists
    merged_df["Position"] = merged_df["IGH_Pos"].fillna(merged_df["Position"]).astype(int)
    merged_df["Allele1"] = merged_df["IGH_Allele"].fillna(merged_df["Allele1"])
    merged_df["Allele2"] = merged_df["IGH_Allele"].fillna(merged_df["Allele2"])

    # Save the updated .bim file
    merged_df[["Chromosome", "SNP", "cM", "Position", "Allele1", "Allele2"]].to_csv(output_bim_file, sep="\t", index=False, header=False)
    print(f"Converted .bim file saved to: {output_bim_file}")


def generate_binary_files(prefix, converted_bim, output_prefix, plink_path="plink"):
    """
    Generate .bim, .bed, .fam files using PLINK.
    
    Arguments:
    - prefix: Prefix of the original files (without extension).
    - converted_bim: Path to the converted .bim file.
    - output_prefix: Path to save the output files.
    - plink_path: Path to plink executable. Default is "plink". Change to "plink2" if needed.
    """
    # Check if the user is using plink2 or plink
    if "plink2" in plink_path.lower():
        cmd = f'{plink_path} --bfile {prefix} --make-bed --out {output_prefix} --bim {converted_bim}'
    else:
        cmd = f'{plink_path} --bfile {prefix} --bim {converted_bim} --make-bed --out {output_prefix}'
    
    print(f"Running command: {cmd}")
    os.system(cmd)
    print(f"Binary files saved to: {output_prefix}")
