import gzip

# ======== Update these paths ========
vcf_in = "GoalPost_hg19_QC_fixed_chr14_region.vcf.gz"   # or GoalPost_hg38_QC_fixed_chr14_region.vcf.gz
mapping_file = "igh_to_hg19.txt"                         # or igh_to_hg38.txt
vcf_out = vcf_in.replace(".vcf.gz", "_custom.vcf")
tsv_out = vcf_in.replace(".vcf.gz", "_custom.tsv")
log_out = vcf_in.replace(".vcf.gz", "_mismatch.log")
# ====================================

# Load mapping into dictionary: {(chrom, pos): (custom_pos, ref_allele)}
pos_map = {}
with open(mapping_file) as f:
    for line in f:
        chrom, pos, ref_old, custom_pos, ref_new = line.strip().split()
        if chrom == "chr14":
            pos_map[("14", pos)] = (custom_pos, ref_new)

# Process VCF
with gzip.open(vcf_in, "rt") as fin, \
     open(vcf_out, "w") as fout_vcf, \
     open(tsv_out, "w") as fout_tsv, \
     open(log_out, "w") as flog:

    for line in fin:
        if line.startswith("##"):
            fout_vcf.write(line)
            continue
        elif line.startswith("#CHROM"):
            header = line.strip().split()
            sample_ids = header[9:]
            fout_vcf.write(line)
            fout_tsv.write("CHR\tCUSTOM_POS\tREF\tALT\t" + "\t".join(sample_ids) + "\n")
            continue

        cols = line.strip().split("\t")
        chrom, pos, id_, ref, alt, qual, filt, info, fmt = cols[:9]
        genotypes = cols[9:]

        key = (chrom, pos)
        if key not in pos_map:
            flog.write(f"SKIP: {chrom}:{pos} not found in mapping\n")
            continue

        new_pos, ref_map = pos_map[key]

        # Match or flip
        if ref == ref_map:
            new_line = [chrom, new_pos, id_, ref, alt, qual, filt, info, fmt] + genotypes
            fout_vcf.write("\t".join(new_line) + "\n")
            fout_tsv.write(f"{chrom}\t{new_pos}\t{ref}\t{alt}\t" + "\t".join(genotypes) + "\n")
        elif alt == ref_map:
            # Flip REF/ALT and genotypes
            flipped_gt = []
            for g in genotypes:
                g = g.replace("0/0", "1/1").replace("1/1", "0/0").replace("0/1", "1/0").replace("1/0", "0/1")
                flipped_gt.append(g)
            new_line = [chrom, new_pos, id_, alt, ref, qual, filt, info, fmt] + flipped_gt
            fout_vcf.write("\t".join(new_line) + "\n")
            fout_tsv.write(f"{chrom}\t{new_pos}\t{alt}\t{ref}\t" + "\t".join(flipped_gt) + "\n")
        else:
            flog.write(f"MISMATCH: {chrom}:{pos} REF={ref} ALT={alt} ≠ expected {ref_map}\n")
