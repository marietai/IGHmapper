# create_shared_position_ranges.py

goalpost_pvar = r"C:\Users\marietai\Desktop\GoalPost\PLINK_QC\LDpruned_data_pgen.pvar"
thousandg_pvar = r"C:\Users\marietai\Desktop\GoalPost\1000G\all_hg38.pvar"
output_file = r"C:\Users\marietai\Desktop\GoalPost\PLINK_QC\shared_positions.txt"

def read_positions(pvar_path):
    print(f"📥 Reading positions from {pvar_path}...")
    positions = set()
    with open(pvar_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue  # skip headers/comments
            fields = line.strip().split()
            if len(fields) > 2 and fields[1].isdigit():
                chrom = fields[0].replace('chr', '')
                pos = int(fields[1])
                positions.add((chrom, pos))
    print(f"✅ Found {len(positions):,} positions.")
    return positions

# Read positions from both datasets
goalpost_positions = read_positions(goalpost_pvar)
thousandg_positions = read_positions(thousandg_pvar)

# Intersect positions
print("🔗 Finding shared positions...")
shared_positions = goalpost_positions.intersection(thousandg_positions)
print(f"✅ Shared positions: {len(shared_positions):,}")

# Write to extract-range file
print(f"💾 Writing output to {output_file}...")
with open(output_file, 'w') as f:
    for chrom, pos in sorted(shared_positions, key=lambda x: (str(x[0]), int(x[1]))):
        f.write(f"{chrom}\t{pos}\t{pos + 1}\n")

print("✅ Done! File is ready for use with --extract range.")
