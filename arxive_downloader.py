# Script to download the newest arXiv source tarballs for a given category.
import os, time, sys
import arxiv

# --- config ---
CATEGORY     = "astro-ph.CO"            # change if you want: e.g. "astro-ph.CO"
NUM_PAPERS   = 200
OUT_DIR      = os.path.join(os.path.dirname(__file__), "Data", f"arx_sources/{CATEGORY[:-3]}/{CATEGORY[-2:]}/sources")
DELAY_SECONDS = 3                       # arXiv politeness (keep it ≥3)

os.makedirs(OUT_DIR, exist_ok=True)

# Build a search: newest first, cap at 200 results so no deep pagination
search = arxiv.Search(
    query=f"cat:{CATEGORY}",
    max_results=NUM_PAPERS,
    sort_by=arxiv.SortCriterion.SubmittedDate,
    sort_order=arxiv.SortOrder.Descending,
)

client = arxiv.Client(delay_seconds=DELAY_SECONDS)

downloaded = 0
seen = 0

while True:
    try:
        for r in client.results(search):
            seen += 1
            arx_id = r.entry_id.rsplit("/", 1)[-1]
            out_path = os.path.join(OUT_DIR, f"{arx_id}.tar.gz")
            if os.path.exists(out_path):
                continue
            try:
                r.download_source(dirpath=OUT_DIR)
                downloaded += 1
                print(f"[{downloaded}/{NUM_PAPERS}] saved {arx_id}")
            except Exception as e:
                print(f"[skip] {arx_id}: {e}", file=sys.stderr)
        break  # finished normally
    except arxiv.UnexpectedEmptyPageError as e:
        # Rare hiccup; quick retry
        print("[warn] empty page from API; retrying in 5s…", file=sys.stderr)
        time.sleep(5)
        continue

print(f"[done] seen={seen}, newly downloaded={downloaded}, saved to {OUT_DIR}")
