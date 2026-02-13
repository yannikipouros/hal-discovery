import pandas as pd
from pathlib import Path
from ..paths import RunContext


def run(ctx: RunContext) -> Path:
    """
    Part 1:
    - Load the InterPro/PFAM family space table from inputs/
    - Keep only rows marked include == 'T'
    - Write a cleaned CSV to outputs/<RUN_ID>/1_protein_families/
    Returns the path to the cleaned CSV.
    """
    in_tsv = ctx.input_dir / "entry-matching-CL0029.tsv"
    if not in_tsv.exists():
        raise FileNotFoundError(f"Missing required input: {in_tsv}")

    df = pd.read_csv(in_tsv, sep="\t")

    # Light sanity check (won't crash if extra columns exist)
    if "Accession" not in df.columns:
        raise ValueError(f"'Accession' column not found in {in_tsv.name}. Columns: {list(df.columns)}")

    if "include" in df.columns:
        df_keep = df.loc[df["include"].astype(str).str.strip() == "T"].copy()
    else:
        # If include column missing, keep everything (but warn via print)
        df_keep = df.copy()
        print("Warning: 'include' column not found; keeping all rows.")

    out_csv = ctx.dirs["1_protein_families"] / "families_included.csv"
    df_keep.to_csv(out_csv, index=False)

    print(f"Loaded families table: {in_tsv.name} ({len(df)} rows)")
    print(f"Kept families: {len(df_keep)}")
    print(f"Wrote: {out_csv}")

    return out_csv
