#!/usr/bin/env python3
import pandas as pd
import glob
import numpy as np
import os
import argparse
from scipy.stats import hmean

# === Argument parsing ===
parser = argparse.ArgumentParser(description="Compute PHBR scores from NetMHC outputs (correct _twice handling, no triple counting).")
parser.add_argument("mhc1_dir", type=str, help="Directory containing MHC-I NetMHCpan result files")
parser.add_argument("mhc2_dir", type=str, help="Directory containing MHC-II NetMHCIIpan result files")
parser.add_argument("--sample_id", type=str, required=True, help="Sample ID for naming output files")
args = parser.parse_args()

mhc1_dir = args.mhc1_dir
mhc2_dir = args.mhc2_dir
sample_id = args.sample_id

# === Function to process MHC files ===
def process_files(file_pattern, label="mhc1"):
    files = glob.glob(file_pattern)
    data = {}

    print(f"🔍 Processing {len(files)} NetMHC result files for {label.upper()}...")

    for f in files:
        try:
            filename = os.path.basename(f).replace(".xls", "")
            parts = filename.split("_")

            # Correct parsing of allele and _twice
            if parts[-1] == "twice":
                twice = True
                allele = "_".join(parts[2:-1])
            else:
                twice = False
                allele = "_".join(parts[2:])

            allele = allele.replace(":", "")  # remove colons safely
            sample_id_local = parts[1]

            print(f"📄 Reading {filename}.xls ➔ Allele: {allele} (Twice: {twice})")

            # Load file
            df = pd.read_csv(f, sep='\t', comment='#', header=1)
            if not {'Pos', 'Peptide'}.issubset(df.columns):
                print(f"⚠️ Skipping {f}: Missing required columns.")
                continue

            df = df[(df['Pos'].between(0, 10)) & (df['Pos'] + df['Peptide'].str.len() > 10)]
            if df.empty:
                print(f"⚠️ Skipping {f}: No valid peptides after filtering.")
                continue

            possible_rank_cols = [col for col in df.columns if 'rank' in col.lower()]
            rank_col = next((col for col in possible_rank_cols if 'el' in col.lower() or '%' in col), df.columns[-3])
            id_col = df.columns[2]

            df = df[[id_col, rank_col]].copy()
            df.columns = ['ID', 'Rank']

            df['ID'] = df['ID'].str.replace('p_', 'p.', regex=False)
            df['gene'] = df['ID'].str.extract(r'^(.*?)_p')[0]
            df['mutation'] = df['ID'].str.extract(r'(p\.[A-Z][0-9]+[A-Z*]+)')[0]

            df = df.dropna(subset=['gene', 'mutation', 'Rank'])
            df['Rank'] = pd.to_numeric(df['Rank'], errors='coerce')
            df = df.dropna(subset=['Rank'])

            if df.empty:
                continue

            grouped = df.groupby(['gene', 'mutation'])['Rank'].min()

            if sample_id_local not in data:
                data[sample_id_local] = {}

            data[sample_id_local][allele + ("_twice" if twice else "")] = grouped

        except Exception as e:
            print(f"❌ Error parsing {f}: {e}")
            continue

    # 🚨 Fix: Remove normal allele if _twice exists
    for sid in data:
        alleles = list(data[sid].keys())
        for allele in alleles:
            if allele.endswith("_twice"):
                normal = allele.replace("_twice", "")
                if normal in data[sid]:
                    print(f"⚠️ Skipping normal {normal} because {allele} exists for {sid}")
                    del data[sid][normal]

    results = []
    for sid, allele_dict in data.items():
        df = pd.DataFrame(allele_dict)
        df = df.reset_index()

        # Expand values properly
        def expanded_ranks(row):
            values = []
            for colname, value in row.items():
                if colname in ['gene', 'mutation']:
                    continue
                if pd.isna(value) or value <= 0:
                    continue
                if colname.endswith("_twice"):
                    values.append(value)
                    values.append(value)  # Count twice
                else:
                    values.append(value)
            return values

        expanded = df.apply(expanded_ranks, axis=1)

        # Calculate PHBR and alleles_used
        df[f'PHBR_{label}'] = expanded.apply(lambda x: hmean(x) if len(x) >= 2 else np.nan)
        df['alleles_used'] = expanded.apply(len)

        # Print expanded values and PHBR
        for idx, values in enumerate(expanded):
            gene = df.loc[idx, 'gene']
            mutation = df.loc[idx, 'mutation']
            print(f"🔬 {sid} ➔ {gene} {mutation}")
            print(f"   ➔ Expanded ranks: {values}")
            print(f"   ➔ Alleles used: {len(values)}")
            if len(values) >= 2:
                print(f"   ➔ PHBR ({label.upper()}): {hmean(values):.3f}")
            else:
                print(f"⚠️ Not enough alleles for PHBR.")

        df[f'PHBR_{label}'] = df[f'PHBR_{label}']
        df['sample_id'] = sid

        results.append(df[['sample_id', 'gene', 'mutation', f'PHBR_{label}', 'alleles_used']])

    if results:
        return pd.concat(results)
    else:
        return pd.DataFrame(columns=['sample_id', 'gene', 'mutation', f'PHBR_{label}', 'alleles_used'])

# === Main ===

# Process MHC-I and MHC-II
mhc1_phbr = process_files(os.path.join(mhc1_dir, "mhc1_*.xls"), label="mhc1")
mhc2_phbr = process_files(os.path.join(mhc2_dir, "mhc2_*.xls"), label="mhc2")

# Save results
if not mhc1_phbr.empty:
    mhc1_phbr.to_csv(f"phbr_{sample_id}_mhc1.csv", index=False)
    print(f"✅ Saved: phbr_{sample_id}_mhc1.csv")

if not mhc2_phbr.empty:
    mhc2_phbr.to_csv(f"phbr_{sample_id}_mhc2.csv", index=False)
    print(f"✅ Saved: phbr_{sample_id}_mhc2.csv")

# Merge if both exist
if not mhc1_phbr.empty and not mhc2_phbr.empty:
    merged = pd.merge(
        mhc1_phbr,
        mhc2_phbr,
        on=['sample_id', 'gene', 'mutation'],
        how='outer',
        suffixes=('_mhc1', '_mhc2')
    )
    merged.to_csv(f"phbr_{sample_id}_merged.csv", index=False)
    print(f"✅ Saved: phbr_{sample_id}_merged.csv")

print("✅ PHBR computation completed successfully.")
