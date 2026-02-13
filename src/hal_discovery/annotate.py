import re
import json
import datetime as dt
from pathlib import Path
from typing import List, Dict, Any

import pandas as pd


# ---------------- FASTA parsing ----------------

_UNIPROT_HEADER_RE = re.compile(
    r"^>(?:sp|tr)\|(?P<acc>[^|]+)\|.*?\s(?P<desc>.*?)\sOS=(?P<os>.+?)\sOX=(?P<ox>\d+)(?:\sGN=(?P<gn>\S+))?.*?$"
)

def parse_uniprot_idmapping_fasta(fasta_path: Path) -> pd.DataFrame:
    """
    Parse UniProt idmapping.fasta output into a DataFrame with:
      structure (accession), function_annotation, organism_scientific,
      organism_taxID, gene_ID, aa_seq
    """
    records: List[Dict[str, Any]] = []

    entry: Dict[str, Any] = {}
    seq_lines: List[str] = []

    with fasta_path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # flush previous
                if entry:
                    entry["aa_seq"] = "".join(seq_lines)
                    records.append(entry)
                    entry = {}
                    seq_lines = []

                m = _UNIPROT_HEADER_RE.match(line)
                if m:
                    entry = {
                        "structure": m.group("acc"),
                        "function_annotation": m.group("desc"),
                        "organism_scientific": m.group("os"),
                        "organism_taxID": m.group("ox"),
                        "gene_ID": m.group("gn") if m.group("gn") else None,
                    }
                else:
                    # fallback: store accession if possible, keep others blank
                    # Example headers can vary; we avoid failing hard.
                    entry = {
                        "structure": None,
                        "function_annotation": None,
                        "organism_scientific": None,
                        "organism_taxID": None,
                        "gene_ID": None,
                    }
            else:
                seq_lines.append(line)

        # flush last
        if entry:
            entry["aa_seq"] = "".join(seq_lines)
            records.append(entry)

    df = pd.DataFrame(records)
    # drop entries without accession
    df = df[df["structure"].notna()].copy()
    return df


# ---------------- filtering ----------------

def build_exclude_pattern(keywords: List[str]) -> re.Pattern:
    # Escape keywords to prevent regex issues, join as OR, case-insensitive
    escaped = [re.escape(k) for k in keywords if k and str(k).strip()]
    if not escaped:
        return re.compile(r"(?!x)x")  # matches nothing
    return re.compile("|".join(escaped), flags=re.IGNORECASE)


def run_annotate(args) -> None:
    hits_csv = Path(args.hits_csv).expanduser().resolve()
    fasta_path = Path(args.idmapping_fasta).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    ts = dt.datetime.now().strftime("%Y%m%d_%H%M%S")

    print(f"Loading hits CSV: {hits_csv.name}")
    df_hits = pd.read_csv(hits_csv)

    if "structure" not in df_hits.columns:
        raise ValueError(f"'structure' column missing from hits CSV: {hits_csv}")

    print(f"Parsing UniProt FASTA: {fasta_path.name}")
    df_uni = parse_uniprot_idmapping_fasta(fasta_path)

    print(f"UniProt entries parsed: {len(df_uni)}")

    # merge
    df_merged = pd.merge(df_hits, df_uni, on="structure", how="left")

    # keyword filter
    pattern = build_exclude_pattern(args.exclude_keywords)
    func = df_merged["function_annotation"].fillna("")

    df_filtered = df_merged[~func.str.contains(pattern)].copy()

    # write outputs
    out_filtered = outdir / f"putativeHals_annotated_and_filtered_{ts}.csv"
    out_ids = outdir / f"putativeHals_annotated_and_filtered_IDsONLY_{ts}.csv"

    df_filtered.to_csv(out_filtered, index=False)
    df_filtered[["structure"]].drop_duplicates().to_csv(out_ids, index=False)

    summary = {
        "hits_csv": str(hits_csv),
        "idmapping_fasta": str(fasta_path),
        "out_filtered_csv": str(out_filtered),
        "out_ids_csv": str(out_ids),
        "n_hits_in": int(len(df_hits)),
        "n_hits_annotated": int(df_merged["function_annotation"].notna().sum()),
        "n_retained_after_filter": int(len(df_filtered)),
        "n_unique_accessions_retained": int(df_filtered["structure"].nunique()),
        "exclude_keywords": list(args.exclude_keywords),
    }

    (outdir / f"summary_annotate_{ts}.json").write_text(json.dumps(summary, indent=2))

    print("\n✅ Annotate complete.")
    print(summary)
