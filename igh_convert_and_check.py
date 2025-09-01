#!/usr/bin/env python3
"""
IGH custom-reference utilities (GRCh38 -> IGH) in one CLI.

Subcommands:
  1) identify-locus:  Read `igh_to_hg38.txt` and print chr/min/max GRCh38 positions covered.
  2) extract-plink:   Use plink2 to extract chr14:from-bp–to-bp to a new prefix.
  3) convert-bim:     Convert a .bim from GRCh38 coords to IGH coords using a map file.
  4) regen-bedfam:    Use plink2 to rebuild .bed/.fam with a custom .bim.
  5) confirm-mapping: Compare a custom .bim against the map and write success/fail lists.

Typical usage (Windows-friendly):
  python igh_convert_and_check.py identify-locus --map C:\path\igh_to_hg38.txt

  plink2 --bfile C:\path\GoalPost_final_QCpassed ^
         --chr 14 --from-bp 105860505 --to-bp 106883718 ^
         --make-bed --out C:\path\GoalPost_IGH_loci

  python igh_convert_and_check.py convert-bim ^
    --bim C:\path\GoalPost_IGH_loci.bim ^
    --map C:\path\igh_to_hg38.txt ^
    --out-bim C:\path\GoalPost_IGH_reference.bim --overwrite-alleles

  plink2 --bfile C:\path\GoalPost_IGH_loci ^
         --bim C:\path\GoalPost_IGH_reference.bim ^
         --make-bed --out C:\path\GoalPost_IGH_reference

  python igh_convert_and_check.py confirm-mapping ^
    --custom-bim C:\path\GoalPost_IGH_reference.bim ^
    --map C:\path\igh_to_hg38.txt ^
    --out-success C:\path\successfully_mapped_variants.txt ^
    --out-fail C:\path\unsuccessfully_mapped_variants.txt

You can also call plink2 from this script with `extract-plink` and `regen-bedfam`
if you pass `--run` (otherwise it just prints the exact command).
"""

from __future__ import annotations
import argparse
import subprocess
from pathlib import Path
import sys
import pandas as pd


# ----------------------------
# utilities
# ----------------------------

def _norm_chr(col: pd.Series) -> pd.Series:
    """Strip leading 'chr' if present; keep as string."""
    return col.astype(str).str.replace("chr", "", regex=False)


def _require_exists(path: Path, kind: str = "file"):
    if not path.exists():
        raise FileNotFoundError(f"{kind.capitalize()} not found: {path}")


def _read_map(map_path: Path) -> pd.DataFrame:
    # Map file columns: Chromosome, GRCh38_Pos, GRCh38_Allele, IGH_Pos, IGH_Allele
    df = pd.read_csv(
        map_path, sep="\t", header=None,
        names=["Chromosome", "GRCh38_Pos", "GRCh38_Allele", "IGH_Pos", "IGH_Allele"]
    )
    df["Chromosome"] = _norm_chr(df["Chromosome"])
    # ensure numeric where appropriate
    df["GRCh38_Pos"] = pd.to_numeric(df["GRCh38_Pos"], errors="coerce").astype("Int64")
    df["IGH_Pos"]     = pd.to_numeric(df["IGH_Pos"], errors="coerce").astype("Int64")
    return df


def _read_bim(bim_path: Path) -> pd.DataFrame:
    # Standard BIM columns
    df = pd.read_csv(
        bim_path, sep="\t", header=None,
        names=["Chromosome", "SNP", "cM", "Position", "Allele1", "Allele2"],
        dtype={"Chromosome": str}
    )
    # enforce string chrom (no "1"->1 surprises)
    df["Chromosome"] = df["Chromosome"].astype(str)
    df["Position"]   = pd.to_numeric(df["Position"], errors="coerce").astype("Int64")
    return df


def _write_bim(df: pd.DataFrame, out_path: Path):
    df.to_csv(out_path, sep="\t", index=False, header=False)


# ----------------------------
# subcommands
# ----------------------------

def cmd_identify_locus(args: argparse.Namespace) -> int:
    map_df = _read_map(Path(args.map))
    if map_df.empty:
        print("Map file is empty.", file=sys.stderr)
        return 1
    chromosome = map_df["Chromosome"].iloc[0]
    # sanity: ensure all in same chromosome
    if (map_df["Chromosome"] != chromosome).any():
        print("Warning: map file spans multiple chromosomes; using overall min/max.", file=sys.stderr)

    min_pos = int(map_df["GRCh38_Pos"].min())
    max_pos = int(map_df["GRCh38_Pos"].max())
    print(f"Loci covered by map: Chromosome {chromosome} from {min_pos} to {max_pos}")
    return 0


def cmd_extract_plink(args: argparse.Namespace) -> int:
    """
    Construct and optionally run plink2 --chr 14 --from-bp ... --to-bp ...
    """
    bfile = Path(args.bfile_prefix)
    _require_exists(Path(str(bfile) + ".bed"))
    _require_exists(Path(str(bfile) + ".bim"))
    _require_exists(Path(str(bfile) + ".fam"))

    out_prefix = Path(args.out_prefix)
    cmd = [
        args.plink2, "--bfile", str(bfile),
        "--chr", str(args.chrom),
        "--from-bp", str(args.from_bp),
        "--to-bp", str(args.to_bp),
        "--make-bed",
        "--out", str(out_prefix),
    ]
    print("PLINK2 command:")
    print(" ".join(cmd))
    if args.run:
        print("Running plink2...")
        completed = subprocess.run(cmd)
        return completed.returncode
    return 0


def cmd_convert_bim(args: argparse.Namespace) -> int:
    """
    Convert a GRCh38 .bim to IGH coordinates using the map.
    - Position -> IGH_Pos (when available)
    - Optionally overwrite Allele1/Allele2 with IGH_Allele (if --overwrite-alleles)
    """
    bim_df = _read_bim(Path(args.bim))
    map_df = _read_map(Path(args.map))

    # If user specified --chrom, filter map to that chrom; else infer from BIM (majority value)
    if args.chrom:
        map_df = map_df[map_df["Chromosome"] == str(args.chrom)]
        if map_df.empty:
            print(f"No rows in map for chrom {args.chrom}", file=sys.stderr)
            return 1
        # also normalize BIM chromosome to that chrom if needed
        bim_df["Chromosome"] = str(args.chrom)
    else:
        # ensure BIM chromosome consistent (many PLINK .bim are "14")
        major_chrom = bim_df["Chromosome"].mode().iat[0]
        map_df = map_df[map_df["Chromosome"] == major_chrom]

    merged = pd.merge(
        bim_df, map_df,
        how="left",
        left_on=["Chromosome", "Position"],
        right_on=["Chromosome", "GRCh38_Pos"]
    )

    if merged.empty:
        print("Error: No matching positions found between .bim and map file.", file=sys.stderr)
        return 1

    # Update position
    merged["Position"] = merged["IGH_Pos"].fillna(merged["Position"]).astype("Int64")

    # Optionally update alleles from map (set both alleles to IGH_Allele like your original snippet)
    if args.overwrite_alleles:
        merged["Allele1"] = merged["IGH_Allele"].fillna(merged["Allele1"])
        merged["Allele2"] = merged["IGH_Allele"].fillna(merged["Allele2"])

    out_cols = ["Chromosome", "SNP", "cM", "Position", "Allele1", "Allele2"]
    out_bim = Path(args.out_bim)
    out_bim.parent.mkdir(parents=True, exist_ok=True)
    _write_bim(merged[out_cols], out_bim)
    print(f"✅ Wrote converted BIM: {out_bim}")

    # Report mapping stats
    n_total = len(bim_df)
    n_mapped = merged["IGH_Pos"].notna().sum()
    print(f"Mapping stats: mapped {n_mapped} / {n_total} ({n_mapped/n_total:.1%})")
    return 0


def cmd_regen_bedfam(args: argparse.Namespace) -> int:
    """
    Rebuild .bed/.fam using an existing .bim override via plink2.
    --bfile-prefix is the original extracted-prefix (with .bed/.fam)
    --bim is the converted .bim
    """
    bfile = Path(args.bfile_prefix)
    _require_exists(Path(str(bfile) + ".bed"))
    _require_exists(Path(str(bfile) + ".fam"))
    _require_exists(Path(args.bim))

    out_prefix = Path(args.out_prefix)
    cmd = [
        args.plink2, "--bfile", str(bfile),
        "--bim", str(Path(args.bim)),
        "--make-bed",
        "--out", str(out_prefix),
    ]
    print("PLINK2 command:")
    print(" ".join(cmd))
    if args.run:
        print("Running plink2...")
        completed = subprocess.run(cmd)
        return completed.returncode
    return 0


def cmd_confirm_mapping(args: argparse.Namespace) -> int:
    """
    Compare a custom .bim (IGH coords) against the map (IGH_Pos) and
    write successfully/unsuccessfully mapped variants.
    """
    custom_bim = _read_bim(Path(args.custom_bim))
    map_df = _read_map(Path(args.map))

    # Use IGH positions for matching
    # If user passed chrom, filter map
    if args.chrom:
        map_df = map_df[map_df["Chromosome"] == str(args.chrom)]

    # Left join on Position == IGH_Pos
    merged = pd.merge(
        custom_bim, map_df[["IGH_Pos"]],
        left_on="Position", right_on="IGH_Pos", how="left"
    )
    success_mask = merged["IGH_Pos"].notna()

    success_df = custom_bim[success_mask]
    fail_df    = custom_bim[~success_mask]

    out_success = Path(args.out_success) if args.out_success else Path("successfully_mapped_variants.txt")
    out_fail    = Path(args.out_fail) if args.out_fail else Path("unsuccessfully_mapped_variants.txt")

    success_df.to_csv(out_success, sep="\t", index=False)
    fail_df.to_csv(out_fail, sep="\t", index=False)

    print(f"✅ Successfully mapped variants: {len(success_df)} -> {out_success}")
    print(f"❗ Unsuccessfully mapped variants: {len(fail_df)} -> {out_fail}")
    return 0


# ----------------------------
# main / argparse
# ----------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="IGH custom-reference utilities (GRCh38 -> IGH)."
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    # identify-locus
    sp = sub.add_parser("identify-locus", help="Print chr/min/max GRCh38 positions covered by the map.")
    sp.add_argument("--map", required=True, help="Path to igh_to_hg38.txt")
    sp.set_defaults(func=cmd_identify_locus)

    # extract-plink
    sp = sub.add_parser("extract-plink", help="Use plink2 to extract a chr region to a new dataset.")
    sp.add_argument("--plink2", default="plink2", help="Path to plink2 executable (default: plink2)")
    sp.add_argument("--bfile-prefix", required=True, help="Input bfile prefix (without extension)")
    sp.add_argument("--chrom", default="14", help="Chromosome to extract (default: 14)")
    sp.add_argument("--from-bp", type=int, required=True, help="Start bp (GRCh38)")
    sp.add_argument("--to-bp",   type=int, required=True, help="End bp (GRCh38)")
    sp.add_argument("--out-prefix", required=True, help="Output bfile prefix")
    sp.add_argument("--run", action="store_true", help="Actually run plink2 (otherwise just print the command)")
    sp.set_defaults(func=cmd_extract_plink)

    # convert-bim
    sp = sub.add_parser("convert-bim", help="Convert a .bim from GRCh38 to IGH coords using the map.")
    sp.add_argument("--bim", required=True, help="Input .bim path (GRCh38 coords)")
    sp.add_argument("--map", required=True, help="Path to igh_to_hg38.txt")
    sp.add_argument("--out-bim", required=True, help="Output .bim path (IGH coords)")
    sp.add_argument("--chrom", default=None, help="Force chromosome label (e.g., 14). If omitted, infer from BIM.")
    sp.add_argument("--overwrite-alleles", action="store_true",
                    help="Overwrite Allele1/Allele2 with IGH_Allele from the map (like your original snippet).")
    sp.set_defaults(func=cmd_convert_bim)

    # regen-bedfam
    sp = sub.add_parser("regen-bedfam", help="Rebuild .bed/.fam with a custom .bim using plink2.")
    sp.add_argument("--plink2", default="plink2", help="Path to plink2 executable (default: plink2)")
    sp.add_argument("--bfile-prefix", required=True, help="Original extracted bfile prefix (has .bed/.fam)")
    sp.add_argument("--bim", required=True, help="Converted .bim path (IGH coords)")
    sp.add_argument("--out-prefix", required=True, help="Output bfile prefix")
    sp.add_argument("--run", action="store_true", help="Actually run plink2 (otherwise just print the command)")
    sp.set_defaults(func=cmd_regen_bedfam)

    # confirm-mapping
    sp = sub.add_parser("confirm-mapping", help="Compare custom .bim against the map; write success/fail lists.")
    sp.add_argument("--custom-bim", required=True, help="Custom .bim in IGH coords")
    sp.add_argument("--map", required=True, help="Path to igh_to_hg38.txt")
    sp.add_argument("--chrom", default=None, help="Chromosome filter (e.g., 14)")
    sp.add_argument("--out-success", default=None, help="Output path for successfully mapped variants")
    sp.add_argument("--out-fail", default=None, help="Output path for unsuccessfully mapped variants")
    sp.set_defaults(func=cmd_confirm_mapping)

    return p


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    try:
        return args.func(args)
    except FileNotFoundError as e:
        print(str(e), file=sys.stderr)
        return 2
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 3


if __name__ == "__main__":
    raise SystemExit(main())
