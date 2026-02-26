#!/usr/bin/env python3
"""
Download GTEx tissue gene-level matrices from recount3 (raw files).

Example:
  python download_gtex_recount3.py --out raw-data --tissues brain heart liver pancreas --n 4

Notes:
- This downloads "gene_sums" matrices (GENCODE v26, G026).
- Files can be large. You'll likely want to process them in chunks downstream.
"""

import argparse
import gzip
import os
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests


RECOUNT3_URL_DEFAULT = "http://duffel.rail.bio/recount3"  # default used by recount3 docs
ORGANISM = "human"
PROJECT_HOME = "data_sources/gtex"
ANN_EXT = "G026"  # GENCODE v26 gene annotation in recount3


def mkdirp(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def download_file(url: str, out_path: Path, sleep: float, session: requests.Session) -> bool:
    if out_path.exists() and out_path.stat().st_size > 0:
        return False
    tmp = out_path.with_suffix(out_path.suffix + ".part")
    with session.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    tmp.replace(out_path)
    if sleep > 0:
        time.sleep(sleep)
    return True


def read_gz_text(path: Path, max_bytes: int = 2_000_000) -> str:
    # Read a small prefix to detect delimiter/header safely
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        return f.read(max_bytes)


def parse_project_table(gz_path: Path) -> List[Dict[str, str]]:
    """
    Parse gtex.recount_project.MD.gz into a list of dict rows.
    Format is tab-delimited in recount3 metadata files (typical).
    """
    import csv

    sample = read_gz_text(gz_path)
    # Heuristic: prefer tab if present in header line
    header_line = sample.splitlines()[0] if sample else ""
    delim = "\t" if "\t" in header_line else ","

    rows: List[Dict[str, str]] = []
    with gzip.open(gz_path, "rt", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter=delim)
        for row in reader:
            rows.append({k: (v or "").strip() for k, v in row.items()})
    if not rows:
        raise RuntimeError(f"Could not parse any rows from {gz_path} (delimiter guess={repr(delim)})")
    return rows


def locate_gene_url(recount3_url: str, project: str) -> str:
    """
    Mirror of recount3::locate_url() logic for type='gene' in data_sources/gtex.

    base_url = file.path(recount3_url, organism, project_home, "gene_sums", last2(project), project)
    file     = paste0(basename(project_home), ".gene_sums.", project, ".", ANN_EXT, ".gz")
    """
    last2 = project[-2:]
    return (
        f"{recount3_url}/{ORGANISM}/{PROJECT_HOME}/gene_sums/{last2}/{project}/"
        f"gtex.gene_sums.{project}.{ANN_EXT}.gz"
    )


def locate_annotation_url(recount3_url: str) -> str:
    return f"{recount3_url}/{ORGANISM}/annotations/gene_sums/{ORGANISM}.gene_sums.{ANN_EXT}.gtf.gz"


def score_match(keyword: str, project: str) -> Tuple[int, int]:
    """
    Higher is better.
    - exact match beats substring
    - shorter project names slightly preferred when ties
    """
    k = keyword.upper()
    p = project.upper()
    if p == k:
        return (10_000, -len(project))
    if k in p:
        # prefer word-boundary-ish matches
        boundary_bonus = 200 if re.search(rf"\b{re.escape(k)}\b", p) else 0
        return (1_000 + boundary_bonus, -len(project))
    return (0, -len(project))


def pick_projects(rows: List[Dict[str, str]], tissues: List[str]) -> Dict[str, str]:
    """
    Map each requested tissue keyword -> best matching recount3 GTEx project.
    """
    # Find the column that contains the project name.
    # In recount3 metadata it's typically "project" or "recount_project.project".
    candidate_cols = ["project", "recount_project.project"]
    col = None
    for c in candidate_cols:
        if c in rows[0]:
            col = c
            break
    if col is None:
        # fallback: pick first col containing "project"
        for c in rows[0].keys():
            if "project" in c.lower():
                col = c
                break
    if col is None:
        raise RuntimeError(f"Could not find a project column in metadata. Columns={list(rows[0].keys())[:20]}")

    projects = sorted({r[col] for r in rows if r.get(col)})
    if not projects:
        raise RuntimeError("No projects found in metadata table.")

    out: Dict[str, str] = {}
    for t in tissues:
        best = None
        best_score = (-1, 0)
        for p in projects:
            sc = score_match(t, p)
            if sc > best_score:
                best_score = sc
                best = p
        if best is None or best_score[0] <= 0:
            raise RuntimeError(
                f"No GTEx project matched tissue keyword {t!r}. "
                f"Try a more specific keyword (e.g. 'heart_left', 'brain_cortex')."
            )
        out[t] = best
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=str, default="raw-data", help="Output directory")
    ap.add_argument("--tissues", nargs="+", default=["brain", "heart", "liver", "pancreas"],
                    help="Tissue keywords (matched against recount3 GTEx project names)")
    ap.add_argument("--recount3-url", type=str, default=RECOUNT3_URL_DEFAULT, help="recount3 base URL")
    ap.add_argument("--sleep", type=float, default=1.0, help="Politeness delay between downloads (seconds)")
    ap.add_argument("--with-annotation", action="store_true", help="Also download the gene annotation GTF")
    args = ap.parse_args()

    out_root = Path(args.out)
    mkdirp(out_root)

    session = requests.Session()
    session.headers.update({"User-Agent": "gtex-recount3-downloader/1.0 (contact: you@example.com)"})

    # 1) Download GTEx project list (once)
    proj_meta_url = f"{args.recount3_url}/{ORGANISM}/{PROJECT_HOME}/metadata/gtex.recount_project.MD.gz"
    proj_meta_path = out_root / "gtex.recount_project.MD.gz"
    print(f"[meta] {proj_meta_url}")
    download_file(proj_meta_url, proj_meta_path, args.sleep, session)

    # 2) Parse and choose projects
    rows = parse_project_table(proj_meta_path)
    mapping = pick_projects(rows, args.tissues)

    print("[plan] tissue keyword -> recount3 GTEx project")
    for k, v in mapping.items():
        print(f"  {k:10s} -> {v}")

    # 3) Download gene matrices
    for tissue_kw, project in mapping.items():
        tissue_dir = out_root / project
        mkdirp(tissue_dir)

        url = locate_gene_url(args.recount3_url, project)
        out_path = tissue_dir / os.path.basename(url)

        print(f"[gene] {project}: {url}")
        try:
            did = download_file(url, out_path, args.sleep, session)
            if did:
                print(f"  saved -> {out_path}")
            else:
                print(f"  exists -> {out_path}")
        except requests.HTTPError as e:
            print(f"[error] failed {project} ({tissue_kw}) : {e}", file=sys.stderr)

    # 4) Optional: annotation
    if args.with_annotation:
        ann_url = locate_annotation_url(args.recount3_url)
        ann_path = out_root / os.path.basename(ann_url)
        print(f"[ann] {ann_url}")
        download_file(ann_url, ann_path, args.sleep, session)
        print(f"  saved -> {ann_path}")

    print("[done]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
