from dataclasses import dataclass
from pathlib import Path
import datetime as dt


@dataclass
class RunContext:
    project_root: Path
    input_dir: Path
    outdir: Path
    run_id: str
    dirs: dict


def make_context(project_root: Path, run_id: str | None = None) -> RunContext:
    project_root = project_root.expanduser().resolve()

    input_dir = project_root / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)

    if run_id is None:
        run_id = dt.datetime.now().strftime("%Y%m%d_%H%M%S")

    outdir = project_root / "outputs" / run_id
    outdir.mkdir(parents=True, exist_ok=True)

    dirs = {
        "1_protein_families": outdir / "1_protein_families",
        "2_retrieved_accessionIDs": outdir / "2_retrieved_accessionIDs",
        "3_processed_accessionIDs": outdir / "3_processed_accessionIDs",
        "4_retrieved_AF2modelsche": outdir / "4_retrieved_AF2modelsche",
        "5_find_2His": outdir / "5_find_2His",
        "6_retrieved_halgogenases_from2HisSite": outdir / "6_retrieved_halgogenases_from2HisSite",
        "7_getUNIPROTannotations_IDmapping": outdir / "7_getUNIPROTannotations_IDmapping",
        "8_TableS2": outdir / "8_TableS2",
    }
    for p in dirs.values():
        p.mkdir(parents=True, exist_ok=True)

    return RunContext(
        project_root=project_root,
        input_dir=input_dir,
        outdir=outdir,
        run_id=run_id,
        dirs=dirs,
    )

