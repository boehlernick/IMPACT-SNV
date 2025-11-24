
#!/usr/bin/env python3
import pandas as pd
import os
import sys

"""
Build per-chromosome FAVOR subset CSVs by joining a FAVOR processed CSV with variant keys TSV.
Join keys: Chromosome.String, Position.String, RefVcf.String, AltVcf.String (from FAVOR),
           chrom, pos, ref, alt (from TSV).
Output: one CSV per chromosome, a favor_subset_map.txt mapping chrom->CSV path,
        and missing_variants.log listing any variants from TSV not found in FAVOR.
If FAVOR_database_example.xlsx is available, column order is read from its first sheet; otherwise
the column order from the FAVOR processed CSV is used.

Usage:
  python build_favor_subset.py <favor_processed_csv> <variants_tsv> <output_dir> [favor_example_xlsx]
"""

def main():
    if len(sys.argv) < 4:
        print("Usage: python build_favor_subset.py <favor_processed_csv> <variants_tsv> <output_dir> [favor_example_xlsx]")
        sys.exit(1)

    favor_csv = sys.argv[1]
    variants_tsv = sys.argv[2]
    output_dir = sys.argv[3]
    favor_example_xlsx = sys.argv[4] if len(sys.argv) > 4 else None

    os.makedirs(output_dir, exist_ok=True)

    # Load FAVOR and variants
    favor = pd.read_csv(favor_csv)
    variants = pd.read_csv(variants_tsv, sep='\t', header=None, names=['chrom','pos','ref','alt'])

    # Normalize types (strip 'chr' if present; ensure strings)
    variants['chrom'] = variants['chrom'].astype(str).str.replace('chr','', regex=False)
    variants['pos'] = variants['pos'].astype(str)
    variants['ref'] = variants['ref'].astype(str)
    variants['alt'] = variants['alt'].astype(str)

    # Validate FAVOR join columns
    required_cols = ['Chromosome.String','Position.String','RefVcf.String','AltVcf.String']
    for col in required_cols:
        if col not in favor.columns:
            print(f"Error: FAVOR file is missing expected column '{col}'. "
                  f"First columns available: {list(favor.columns)[:10]} ...")
            sys.exit(2)

    # Normalize FAVOR key types
    favor['Chromosome.String'] = favor['Chromosome.String'].astype(str).str.replace('chr','', regex=False)
    favor['Position.String'] = favor['Position.String'].astype(str)
    favor['RefVcf.String'] = favor['RefVcf.String'].astype(str)
    favor['AltVcf.String'] = favor['AltVcf.String'].astype(str)

    # Perform left merge from variant keys to FAVOR
    merged = variants.merge(
        favor,
        left_on=['chrom','pos','ref','alt'],
        right_on=['Chromosome.String','Position.String','RefVcf.String','AltVcf.String'],
        how='left',
        indicator=True
    )

    # Log missing
    missing = merged[merged['_merge'] != 'both'][['chrom','pos','ref','alt']]
    if len(missing):
        missing_path = os.path.join(output_dir, 'missing_variants.log')
        missing.to_csv(missing_path, index=False)
        print(f"Warning: {len(missing)} variants not found in FAVOR. See {missing_path}")

    # Keep only matches
    found = merged[merged['_merge'] == 'both']

    # Column order: prefer FAVOR example xlsx if provided; else use FAVOR CSV order
    if favor_example_xlsx and os.path.exists(favor_example_xlsx):
        try:
            example = pd.read_excel(favor_example_xlsx, engine='openpyxl')
            header_order = list(example.columns)
            # Ensure all are present; if not, fall back to FAVOR CSV order
            if not all(col in found.columns for col in header_order):
                header_order = list(favor.columns)
                print("Info: Some example columns not present in processed FAVOR; using FAVOR CSV column order.")
        except Exception as e:
            print(f"Info: Could not read example header from '{favor_example_xlsx}' ({e}); using FAVOR CSV column order.")
            header_order = list(favor.columns)
    else:
        header_order = list(favor.columns)

    # Reorder to chosen header order
    found = found[header_order]

    # Group and write one file per chromosome (Chromosome.String)
    for chrom, sub in found.groupby('Chromosome.String', sort=True):
        out_csv = os.path.join(output_dir, f"favor_chr{chrom}.csv")
        sub.to_csv(out_csv, index=False)
        print(f"Wrote {out_csv} ({len(sub)} rows)")

    # Write mapping file
    map_path = os.path.join(output_dir, 'favor_subset_map.txt')
    with open(map_path, 'w') as fh:
        for chrom in sorted(found['Chromosome.String'].unique(), key=lambda x: (len(str(x)), str(x))):
            fh.write(f"{chrom}\t{os.path.join(output_dir, f'favor_chr{chrom}.csv')}\n")
    print(f"Wrote mapping to {map_path}")

if __name__ == '__main__':
    main()
