import json
import datetime as dt
from pathlib import Path
from typing import Dict, Tuple, List, Any

import pandas as pd
from Bio.PDB import PDBParser


# ----------------------- geometry utils -----------------------

def euclidean_distance(a, b) -> float:
    """Euclidean distance between two 3D points (iterables length 3)."""
    return ((a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2) ** 0.5


# ----------------------- HIS pairing (notebook-like) -----------------------

def extract_his_coordinates(structure):
    """
    Extract coordinates for HIS residues:
      - NE2 (side chain)
      - N and O (backbone)
    Keys are (chain_id, resseq) to avoid collisions across chains.
    """
    his_ne2 = {}
    his_o = {}
    his_n = {}

    for model in structure:
        for chain in model:
            chain_id = chain.id
            for residue in chain:
                if residue.get_resname() == "HIS":
                    resseq = residue.get_id()[1]
                    key = (chain_id, resseq)
                    # direct lookup is faster/cleaner than looping atoms
                    if "NE2" in residue:
                        his_ne2[key] = residue["NE2"].get_coord()
                    if "O" in residue:
                        his_o[key] = residue["O"].get_coord()
                    if "N" in residue:
                        his_n[key] = residue["N"].get_coord()

    return his_ne2, his_o, his_n


def find_his_pairs(his_ne2, his_o, his_n, cutoff_ne2=4.0, cutoff_hbond=4.0):
    """
    Find pairs of HIS residues using your original geometric criteria.
    Pairs are returned as list of ((chain, resi), (chain, resi)).
    Requires all needed atoms present for both residues.
    """
    keys = sorted(set(his_ne2.keys()) & set(his_o.keys()) & set(his_n.keys()))
    pairs = []

    for idx_i in range(len(keys)):
        for idx_j in range(idx_i + 1, len(keys)):
            i = keys[idx_i]
            j = keys[idx_j]

            dist1 = euclidean_distance(his_ne2[i], his_ne2[j])
            if dist1 >= cutoff_ne2:
                continue

            dist2 = euclidean_distance(his_o[i], his_n[j])
            dist3 = euclidean_distance(his_n[i], his_o[j])

            if dist2 < cutoff_hbond and dist3 < cutoff_hbond:
                pairs.append((i, j))

    return pairs


# ----------------------- closest residue distances (notebook-like) -----------------------

KEY_ATOMS = {
    "ASN": ["ND2", "OD1"],
    "GLN": ["NE2", "OE1", "OE2"],
    "CYS": ["SG"],
    "HIS": ["NE2", "ND1"],
    "MET": ["SD"],
    "TYR": ["OH"],
    "TRP": ["NE1"],
    "ASP": ["OD1", "OD2", "CB", "CG"],
    "GLU": ["OE1", "OE2", "CB", "CG", "CD"],
    "ALA": ["CB"],
    "GLY": ["CA"],
    "PHE": ["CD1", "CD2", "CG", "CE1", "CE2", "CZ"],
    "ILE": ["CG1", "CG2", "CD1"],
    "LEU": ["CG", "CD1", "CD2"],
    "VAL": ["CG1", "CG2"],
    "LYS": ["CG", "CD", "CE", "NZ"],
    "PRO": ["CB", "CG", "CD"],
    "SER": ["OG"],
    "THR": ["OG1", "CG2"],
    "ARG": ["CZ", "NE", "NH1", "NH2"],
}
RES_TYPES = list(KEY_ATOMS.keys())


def _get_hx_resname(structure, chain_id: str, resseq: int, hx_offset: int) -> str:
    """
    Notebook-like HX: residue at (resseq + hx_offset) in the same chain, if present.
    Returns 'NA' if missing.
    """
    target = (" ", resseq + hx_offset, " ")
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            try:
                return chain[target].get_resname()
            except Exception:
                return "NA"
    return "NA"


def find_closest_residues(structure, his_pairs, his_ne2, protein_length, pdb_path, hx_offset: int):
    """
    For each HIS pair:
      - compute midpoint between NE2 atoms
      - compute closest distances from midpoint to key atoms of candidate residue types
      - determine HX residue (two residues after first HIS) when possible
    Returns a list of row dicts (one per 2-His site).
    """
    structure_id = Path(pdb_path).stem
    rows = []

    for pair in his_pairs:
        (chainA, resiA), (chainB, resiB) = pair

        midpoint = (
            (his_ne2[pair[0]][0] + his_ne2[pair[1]][0]) / 2,
            (his_ne2[pair[0]][1] + his_ne2[pair[1]][1]) / 2,
            (his_ne2[pair[0]][2] + his_ne2[pair[1]][2]) / 2,
        )

        # IMPORTANT: reset per pair (your earlier bug cause)
        closest = {res: float("inf") for res in RES_TYPES}

        # HX residue (optional use later; we still compute it)
        hx_res = _get_hx_resname(structure, chainA, resiA, hx_offset=hx_offset)

        for model in structure:
            for chain in model:
                for residue in chain:
                    resseq = residue.get_id()[1]
                    if (chain.id, resseq) in [pair[0], pair[1]]:
                        continue

                    resname = residue.get_resname()
                    if resname not in KEY_ATOMS:
                        continue

                    for atom_name in KEY_ATOMS[resname]:
                        if atom_name in residue:
                            d = euclidean_distance(midpoint, residue[atom_name].get_coord())
                            if d < closest[resname]:
                                closest[resname] = d

        # convert inf -> blank for writing CSV, but keep numeric for filtering
        row = {
            "structure": structure_id,
            "chainA": chainA,
            "hisA": resiA,
            "chainB": chainB,
            "hisB": resiB,
            "protein_length": protein_length,
            "HXres": hx_res,
        }
        for res in RES_TYPES:
            row[f"closest_{res}"] = None if closest[res] == float("inf") else float(closest[res])

        rows.append(row)

    return rows


def main_function(pdb_path: str, cutoff_ne2: float, cutoff_hbond: float, hx_offset: int) -> List[Dict[str, Any]]:
    """Parse PDB, find His pairs, compute distances and context features."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", pdb_path)

    protein_length = len(list(structure.get_residues()))

    his_ne2, his_o, his_n = extract_his_coordinates(structure)
    his_pairs = find_his_pairs(his_ne2, his_o, his_n, cutoff_ne2=cutoff_ne2, cutoff_hbond=cutoff_hbond)

    if not his_pairs:
        return []

    return find_closest_residues(structure, his_pairs, his_ne2, protein_length, pdb_path, hx_offset=hx_offset)


# ----------------------- filtering (notebook-like but parameterized) -----------------------

def _coalesce_numeric(s):
    return pd.to_numeric(s, errors="coerce")


def run_mine(args) -> None:
    """
    hal-discovery mine
      Input: folder of PDB files (named <UniProt>.pdb)
      Output: 2-His site table + filtered hit tables + summary
    """
    pdb_dir = Path(args.pdb_dir).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    ts = dt.datetime.now().strftime("%Y%m%d_%H%M%S")

    pdb_files = sorted(pdb_dir.glob("*.pdb"))
    if args.max_pdbs is not None:
        pdb_files = pdb_files[: args.max_pdbs]

    print(f"Scanning {len(pdb_files)} PDB files in: {pdb_dir}")

    all_rows: List[Dict[str, Any]] = []
    no_pairs: List[str] = []

    for i, pdb_path in enumerate(pdb_files, start=1):
        if i == 1 or i % args.progress_every == 0 or i == len(pdb_files):
            print(f"[{i}/{len(pdb_files)}]")

        try:
            rows = main_function(
                str(pdb_path),
                cutoff_ne2=args.his_ne2_cutoff,
                cutoff_hbond=args.his_hbond_cutoff,
                hx_offset=args.hx_offset,
            )
        except Exception:
            # if a PDB is malformed, skip it
            continue

        if not rows:
            no_pairs.append(pdb_path.stem)
            continue

        all_rows.extend(rows)

    # Write raw 2-His site table
    df_2his = pd.DataFrame(all_rows)
    out_2his = outdir / f"2his_sites_all_{ts}.csv"
    df_2his.to_csv(out_2his, index=False)

    out_none = outdir / f"no_pairs_{ts}.txt"
    out_none.write_text("\n".join(no_pairs) + ("\n" if no_pairs else ""))

    # If empty, still write summary and exit
    if df_2his.empty:
        summary = {
            "total_structures": len(pdb_files),
            "total_2his_sites": 0,
            "proteins_with_2his_sites": 0,
            "sites_2his_aspglu": 0,
            "sites_putative_halogenase": 0,
            "raw_sites_csv": str(out_2his),
            "no_pairs_txt": str(out_none),
        }
        (outdir / f"summary_{ts}.json").write_text(json.dumps(summary, indent=2))
        print("No 2-His sites found.")
        print(summary)
        return

    # Ensure numeric columns
    for res in RES_TYPES:
        col = f"closest_{res}"
        if col in df_2his.columns:
            df_2his[col] = _coalesce_numeric(df_2his[col])

    df_2his["protein_length"] = _coalesce_numeric(df_2his["protein_length"])

    # ---------- Facial triad filter (your notebook logic, fixed parentheses) ----------
    # default cutoff: <5 Å for ASP/GLU, and length > 200
    df_facial = df_2his[
        (
            (df_2his["closest_ASP"] <= args.asp_glu_cutoff)
            | (df_2his["closest_GLU"] <= args.asp_glu_cutoff)
        )
        & (df_2his["protein_length"] > args.triad_min_length)
    ].copy()


    out_facial = outdir / f"facial_triad_sites_{ts}.csv"
    df_facial.to_csv(out_facial, index=False)

    # ---------- Putative halogenase filter (your notebook logic, parameterized) ----------
    # “no 3rd ligand” by requiring acids farther than acid_exclude (default 5.5)
    # exclude other ligands by requiring > exclude_ligand_dist (default 4)
    # require Ala/Gly close (default < 7)
    # optional HX motif requirement (default ON), HX allowed residues default ALA/GLY

    # start with base mask
    mask = (
        (df_2his["closest_ASP"] > args.acid_exclude_cutoff)
        & (df_2his["closest_GLU"] > args.acid_exclude_cutoff)
        & (df_2his["closest_ASN"] > args.exclude_ligand_cutoff)
        & (df_2his["closest_GLN"] > args.exclude_ligand_cutoff)
        & (df_2his["closest_CYS"] > args.exclude_ligand_cutoff)
        & (df_2his["closest_MET"] > args.exclude_ligand_cutoff)
        & (df_2his["closest_HIS"] > args.exclude_ligand_cutoff)
        & (df_2his["closest_TRP"] > args.exclude_ligand_cutoff)
        & (df_2his["closest_TYR"] > args.exclude_ligand_cutoff)
        & (
            (df_2his["closest_ALA"] < args.ala_gly_presence_cutoff)
            | (df_2his["closest_GLY"] < args.ala_gly_presence_cutoff)
        )
        & (df_2his["protein_length"] > args.min_length)
    )

    if args.require_hx:
        mask = mask & (df_2his["HXres"].isin(args.hx_allowed))

    df_hals = df_2his[mask].copy()

    out_hals = outdir / f"putative_halogenases_{ts}.csv"
    df_hals.to_csv(out_hals, index=False)

    out_ids = outdir / f"putative_halogenase_accessions_{ts}.txt"
    df_hals["structure"].drop_duplicates().to_csv(out_ids, index=False, header=False)

    # ---------- printed + saved summary (your preferred wording) ----------
    n_sites_total = len(df_2his)
    n_proteins_total = df_2his["structure"].nunique()
    n_sites_facial = len(df_facial)
    n_sites_hals = len(df_hals)

    summary = {
        "total_structures": len(pdb_files),
        "total_2his_sites": int(n_sites_total),
        "proteins_with_2his_sites": int(n_proteins_total),
        "sites_2his_aspglu": int(n_sites_facial),
        "sites_putative_halogenase": int(n_sites_hals),
        "raw_sites_csv": str(out_2his),
        "facial_triad_csv": str(out_facial),
        "putative_halogenases_csv": str(out_hals),
        "putative_halogenase_ids": str(out_ids),
        "no_pairs_txt": str(out_none),
    }
    (outdir / f"summary_{ts}.json").write_text(json.dumps(summary, indent=2))

    print(
        f"\nThere are {n_sites_total} total 2-His sites identified, in {n_proteins_total} proteins\n"
        f"of which {n_sites_facial} sites are 2-His-1-Asp/Glu,\n"
        f"and {n_sites_hals} are 2-His sites that lack a third metal-coordinating residue and instead contain a nearby Ala/Gly.\n"
    )
    print("Done.")
    print(summary)
