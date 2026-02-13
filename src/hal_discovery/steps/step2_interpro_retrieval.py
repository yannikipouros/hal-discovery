import json
import time
from pathlib import Path
from urllib import request
from urllib.error import HTTPError

import pandas as pd

from ..paths import RunContext


def _base_url_for_id(fam_id: str) -> str:
    fam_id = fam_id.strip()
    if fam_id.startswith("I"):
        return f"https://www.ebi.ac.uk/interpro/api/protein/UniProt/entry/InterPro/{fam_id}/?page_size=200"
    if fam_id.startswith("P"):
        return f"https://www.ebi.ac.uk/interpro/api/protein/UniProt/entry/PFAM/{fam_id}/?page_size=200"
    if fam_id.startswith("S"):
        return f"https://www.ebi.ac.uk/interpro/api/protein/UniProt/entry/ssf/{fam_id}/?page_size=200"
    raise ValueError(f"Unsupported family ID: {fam_id}")


def retrieve_accessions_for_family(fam_id: str, out_dir: Path, polite_delay: float = 0.5) -> Path:
    """
    Retrieve UniProt accessions for one InterPro/PFAM/SSF family ID via InterPro API.
    Writes: <fam_id>_accessions_retrieved.csv
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{fam_id}_accessions_retrieved.csv"

    # fresh file each run for that family
    out_path.write_text("AccessionID,Name,Length,Gene,in_alphafold\n")

    next_page = _base_url_for_id(fam_id)

    while next_page:
        try:
            req = request.Request(next_page, headers={"Accept": "application/json"})
            res = request.urlopen(req)

            # InterPro sometimes responds 204 when no content
            if getattr(res, "status", None) == 204:
                break

            payload = json.loads(res.read().decode())
            next_page = payload.get("next")

            results = payload.get("results", [])
            if not results:
                break

            with out_path.open("a") as f:
                for item in results:
                    md = item.get("metadata", {})
                    f.write(
                        f"{md.get('accession','N/A')},"
                        f"{md.get('name','N/A')},"
                        f"{md.get('length','N/A')},"
                        f"{md.get('gene','N/A')},"
                        f"{md.get('in_alphafold','N/A')}\n"
                    )

            time.sleep(polite_delay)

        except HTTPError as e:
            # 408 timeout: wait and retry same page
            if getattr(e, "code", None) == 408:
                time.sleep(61)
                continue
            # other HTTP errors -> stop this family
            break

        except json.JSONDecodeError:
            break

    return out_path


def run(ctx: RunContext, families_csv: Path) -> list[Path]:
    """
    Part 2:
    - Read families_included.csv from Part 1
    - For each Accession ID, retrieve accessions from InterPro API
    - Write one CSV per family into outputs/<RUN_ID>/2_retrieved_accessionIDs/
    Returns list of output file paths.
    """
    df = pd.read_csv(families_csv)

    if "Accession" not in df.columns:
        raise ValueError(
            f"'Accession' column missing from {families_csv.name}. Columns: {list(df.columns)}"
        )

    out_dir = ctx.dirs["2_retrieved_accessionIDs"]
    out_files: list[Path] = []

    fam_ids = df["Accession"].astype(str).str.strip().tolist()

    print(f"Retrieving accessions for {len(fam_ids)} families…")

    for idx, fam_id in enumerate(fam_ids, start=1):
        expected_file = out_dir / f"{fam_id}_accessions_retrieved.csv"

        # Skip if already downloaded
        if expected_file.exists() and expected_file.stat().st_size > 0:
            if idx == 1 or idx % 10 == 0 or idx == len(fam_ids):
                print(f"[{idx}/{len(fam_ids)}] SKIP {fam_id} (already exists)")
            out_files.append(expected_file)
            continue

        out_path = retrieve_accessions_for_family(
            fam_id,
            out_dir=out_dir,
            polite_delay=0.5
        )
        out_files.append(out_path)

        if idx == 1 or idx % 10 == 0 or idx == len(fam_ids):
            print(f"[{idx}/{len(fam_ids)}] wrote: {out_path.name}")

    return out_files

