# intersect_by_position.py

import time

goalpost_pvar = r"C:\Users\marietai\Desktop\GoalPost\PLINK_QC\LDpruned_data_pgen.pvar"
thousandg_pvar = r"C:\Users\marietai\Desktop\GoalPost\1000G\all_hg38.pvar"
output_file = r"C:\Users\marietai\Desktop\GoalPost\PLINK_QC\shared_snps_by_pos.txt"

def read_variant_positions(pvar_path, label):
    print(f"🔍 Reading positions from {label}...")
    positions = set()
    with open(pvar_path, 'r') as f:
        next(f)
        for line in f:
            fields = line.strip().split()
            if len(fields) > 3:
                chrom = fields[0].replace('chr', '')  # Normalize
                pos = fields[1]
                positions.add((chrom, pos))
    print(f"✅ Done reading {label}: {len(positions):,} positions found.")
    return positions

def get_ids_by_position(pvar_path, position_set, label):
    print(f"📦 Extracting variant IDs from {label}...")
    ids = set()
    with open(pvar_path, 'r') as f:
        next(f)
        for line in f:
            fields = line.strip().split()
            if len(fields) > 3:
                chrom = fields[0].replace('chr', '')
                pos = fields[1]
                snp_id = fields[2]
                if (chrom, pos) in position_set:
                    ids.add(snp_id)
    print(f"✅ Extracted {len(ids):,} variant IDs from {label}.")
    return ids

start_time = time.time()

# Step 1: Read CHR+POS
goalpost_pos = read_variant_positions(goalpost_pvar, "GoalPost")
thousandg_pos = read_variant_positions(thousandg_pvar, "1000 Genomes")

# Step 2: Intersect positions
print("🔗 Finding shared positions...")
shared_pos = goalpost_pos.intersection(thousandg_pos)
print(f"✅ Shared positions: {len(shared_pos):,}")

# Step 3: Get SNP IDs for shared positions
goalpost_ids = get_ids_by_position(goalpost_pvar, shared_pos, "GoalPost")
thousandg_ids = get_ids_by_position(thousandg_pvar, shared_pos, "1000 Genomes")

# Keep all variant IDs from GoalPost at shared positions
print("🔗 Extracting SNPs from GoalPost at shared positions...")
goalpost_shared_ids = get_ids_by_position(goalpost_pvar, shared_pos, "GoalPost")

# Save these to file for PLINK
print(f"💾 Writing output to {output_file}...")
with open(output_file, 'w') as f:
    for snp_id in sorted(goalpost_shared_ids):
        f.write(snp_id + '\n')

print(f"✅ Done! {len(goalpost_shared_ids)} SNPs written based on shared CHR+POS.")
end_time = time.time()
print(f"\n⏱️ Finished in {round(end_time - start_time, 2)} seconds.")
