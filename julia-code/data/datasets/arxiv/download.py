# Download arXiv source tarballs for categories, constrained to a submittedDate window.
import os, time, sys
import argparse
from datetime import datetime, timezone
import arxiv

def parse_iso_to_utc(dt_str: str) -> datetime:
    """
    Accepts:
      - 'YYYY-MM-DD'
      - 'YYYY-MM-DDTHH:MM'
      - 'YYYY-MM-DDTHH:MM:SS'
    If no timezone is provided, it's interpreted as LOCAL time and converted to UTC.
    """
    dt = datetime.fromisoformat(dt_str)
    if dt.tzinfo is None:
        # local time -> UTC
        dt = dt.astimezone()
    return dt.astimezone(timezone.utc)

def arxiv_utc_fmt(dt: datetime) -> str:
    # arXiv wants YYYYMMDDHHMM in GMT/UTC
    return dt.strftime("%Y%m%d%H%M")

# ---------------- CLI ----------------
parser = argparse.ArgumentParser()
parser.add_argument("--n", type=int, default=100, help="Max number of papers to fetch per category (cap).")
parser.add_argument("--from", dest="dt_from", required=True,
                    help="Window start (ISO). e.g. 2026-01-01 or 2026-01-01T00:00")
parser.add_argument("--to", dest="dt_to", required=True,
                    help="Window end (ISO). e.g. 2026-01-29 or 2026-01-29T23:59")
parser.add_argument("--sort", choices=["submitted", "relevance"], default="submitted",
                    help="Sort within the window.")
args = parser.parse_args()

NUM_PAPERS    = args.n
DELAY_SECONDS = 5
MAX_RETRIES   = 5

dt_from_utc = parse_iso_to_utc(args.dt_from)
dt_to_utc   = parse_iso_to_utc(args.dt_to)
if dt_to_utc <= dt_from_utc:
    raise ValueError("--to must be after --from")

from_str = arxiv_utc_fmt(dt_from_utc)
to_str   = arxiv_utc_fmt(dt_to_utc)

# ---------------- categories ----------------
with open("categories.txt") as f:
    CATEGORIES = [line.strip() for line in f if line.strip() and not line.strip().startswith("#")]

client = arxiv.Client(delay_seconds=DELAY_SECONDS)

for CATEGORY in CATEGORIES:
    print(f"\nCategory [{CATEGORY}] window UTC: [{from_str} TO {to_str}]")
    CAT, SUBCAT = CATEGORY.split(".")
    OUT_DIR = os.path.join(os.path.dirname(__file__), "raw-data", CAT, SUBCAT)
    os.makedirs(OUT_DIR, exist_ok=True)

    # Date-window constraint: submittedDate:[FROM TO TO]
    # (This is the official API syntax; it must live inside the query string.)
    query = f"cat:{CATEGORY} AND submittedDate:[{from_str} TO {to_str}]"

    sort_by = arxiv.SortCriterion.SubmittedDate if args.sort == "submitted" else arxiv.SortCriterion.Relevance

    search = arxiv.Search(
        query=query,
        max_results=NUM_PAPERS,
        sort_by=sort_by,
        sort_order=arxiv.SortOrder.Descending,
    )

    downloaded = 0
    seen = 0
    retries = 0

    while retries < MAX_RETRIES:
        try:
            for r in client.results(search):
                seen += 1
                arxiv_id = r.entry_id.rsplit("/", 1)[-1]
                out_path = os.path.join(OUT_DIR, f"{arxiv_id}.tar.gz")
                if os.path.exists(out_path):
                    continue
                try:
                    r.download_source(dirpath=OUT_DIR)
                    downloaded += 1
                    print(f"[{downloaded}] saved {arxiv_id}", end="\r")
                except Exception as e:
                    print(f"\n[skip] {arxiv_id}: {e}", file=sys.stderr)

            break  # success
        except arxiv.UnexpectedEmptyPageError:
            retries += 1
            print("\n[warn] empty page from API; retrying in 5s…", file=sys.stderr)
            time.sleep(5)

    print(f"\n[done] seen={seen}, newly downloaded={downloaded}, saved to {OUT_DIR}")
