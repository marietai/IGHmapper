import pandas as pd
import os

import pandas as pd

def identify_build(bim_file):
    # Load the .bim file
    bim_df = pd.read_csv(bim_file, sep="\t", header=None, 
                         names=["Chromosome", "SNP", "cM", "Position", "Allele1", "Allele2"],
                         dtype={"Chromosome": str})

    # Reference positions based on known build differences
    grch37_ranges = {
        "1": 249250621, "2": 243199373, "3": 198022430,
        "4": 191154276, "5": 180915260, "6": 171115067,
        "7": 159138663, "8": 146364022, "9": 141213431,
        "10": 135534747, "11": 135006516, "12": 133851895,
        "13": 115169878, "14": 107349540, "15": 102531392,
        "16": 90354753, "17": 81195210, "18": 78077248,
        "19": 59128983, "20": 63025520, "21": 48129895,
        "22": 51304566, "X": 155270560, "Y": 59373566
    }
    
    grch38_ranges = {
        "1": 248956422, "2": 242193529, "3": 198295559,
        "4": 190214555, "5": 181538259, "6": 170805979,
        "7": 159345973, "8": 145138636, "9": 138394717,
        "10": 133797422, "11": 135086622, "12": 133275309,
        "13": 114364328, "14": 107043718, "15": 101991189,
        "16": 90338345, "17": 83257441, "18": 80373285,
        "19": 58617616, "20": 64444167, "21": 46709983,
        "22": 50818468, "X": 156040895, "Y": 57227415
    }

    grch37_score = 0
    grch38_score = 0
    total_chromosomes = 0

    chromosomes = bim_df["Chromosome"].unique()

    for chrom in chromosomes:
        if chrom not in grch37_ranges or chrom not in grch38_ranges:
            continue
        
        total_chromosomes += 1
        chr_df = bim_df[bim_df["Chromosome"] == chrom]
        max_pos = chr_df["Position"].max()

        if max_pos <= grch37_ranges[chrom]:
            grch37_score += 1
        if max_pos <= grch38_ranges[chrom]:
            grch38_score += 1

    if total_chromosomes == 0:
        return "Unknown"

    grch37_ratio = grch37_score / total_chromosomes
    grch38_ratio = grch38_score / total_chromosomes

    if grch38_ratio > grch37_ratio:
        return "GRCh38"
    elif grch37_ratio > grch38_ratio:
        return "GRCh37"
    else:
        return "Unknown"


def extract_region(bim_file, bed_file, fam_file, output_prefix):
    # Extracting relevant region: chr14:105000000-107500000
    command = f"plink2 --bfile {bim_file} --chr 14 --from-bp 105000000 --to-bp 107500000 --make-bed --out {output_prefix}"
    os.system(command)
    print(f"Extracted region saved to: {output_prefix}")
