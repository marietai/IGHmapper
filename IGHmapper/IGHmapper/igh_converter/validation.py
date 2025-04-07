import os
import pandas as pd

def validate_mapping(custom_bim_file, build_type, output_success=None, output_failure=None):
    # Automatically load the reference file based on build type
    base_path = os.path.dirname(__file__)

    if build_type == "GRCh37":
        map_file = os.path.join(base_path, "igh_to_hg19.txt")
    elif build_type == "GRCh38":
        map_file = os.path.join(base_path, "igh_to_hg38.txt")
    else:
        raise ValueError("Invalid build type specified. Use 'GRCh37' or 'GRCh38'.")

    # Load the embedded mapping file
    map_df = pd.read_csv(map_file, sep="\t", header=None, 
                         names=["Chromosome", "GRCh38_Pos", "GRCh38_Allele", "IGH_Pos", "IGH_Allele"])
    
    # Remove "chr" prefix if it exists and ensure data types are correct
    map_df["Chromosome"] = map_df["Chromosome"].str.replace("chr", "", regex=False)
    map_df["IGH_Pos"] = map_df["IGH_Pos"].astype(int)

    # Load the .bim file
    custom_bim_df = pd.read_csv(custom_bim_file, sep="\t", header=None, 
                                names=["Chromosome", "SNP", "cM", "Position", "Allele1", "Allele2"])

    # Ensure Position is integer
    custom_bim_df["Position"] = custom_bim_df["Position"].astype(int)

    # Perform the merge based on positions
    successfully_mapped = pd.merge(custom_bim_df, map_df, left_on=["Position"], right_on=["IGH_Pos"])
    unsuccessfully_mapped = custom_bim_df[~custom_bim_df["Position"].isin(successfully_mapped["Position"])]

    # Save results
    if output_success:
        successfully_mapped.to_csv(output_success, sep="\t", index=False)
    if output_failure:
        unsuccessfully_mapped.to_csv(output_failure, sep="\t", index=False)

    print("Validation complete.")
    print(f"Successfully mapped variants saved to: {output_success}")
    print(f"Unsuccessfully mapped variants saved to: {output_failure}")
