# create_sample_labels.py

fam_file = r"C:\Users\marietai\Desktop\GoalPost\PLINK_QC\Merged_GoalPost_1000G.fam"
output_file = r"C:\Users\marietai\Desktop\GoalPost\PLINK_QC\sample_labels.txt"

with open(fam_file, 'r') as f, open(output_file, 'w') as out:
    out.write("IID\tgroup\n")
    for line in f:
        fields = line.strip().split()
        iid = fields[1]
        group = "1000G" if iid.startswith("NA") else "GoalPost"
        out.write(f"{iid}\t{group}\n")

print(f"✅ Sample label file created: {output_file}")
