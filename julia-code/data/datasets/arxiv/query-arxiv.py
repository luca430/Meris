# Script to download the newest arXiv source tarballs for a given category.
import os, time, sys
import arxiv
import argparse

# Simple parser
parser = argparse.ArgumentParser()
parser.add_argument("--n", type=int, default=10, help="Number of papers to fetch.")
args = parser.parse_args()
# Specify parameters
NUM_PAPERS    = args.n   # no. of papers to fetch
DELAY_SECONDS = 5        # arXiv politeness (keep it ≥3)

# Specify/load categories of interest
f = open("categories.txt")
CATEGORIES = f.read().splitlines()
f.close()

# Define arXiv client
client = arxiv.Client(delay_seconds=DELAY_SECONDS)

# For each of the categories, fetch NUM_PAPERS papers
# Sorting is done by 'relevance' [https://info.arxiv.org/help/api/user-manual.html#sort]
for CATEGORY in CATEGORIES:
    print(f"Getting most RELEVANT papers from [{CATEGORY}]")
    CAT, SUBCAT = CATEGORY.split(".")
    DIRNAME = f"{CAT}/{SUBCAT}/"
    OUT_DIR = os.path.join(os.path.dirname(__file__), "raw-data/", DIRNAME)

    # Make directory if it does not exist
    os.makedirs(OUT_DIR, exist_ok=True)

    # Build a search: most 'relevant' first
    search = arxiv.Search(
        query=f"cat:{CATEGORY}",
        max_results=NUM_PAPERS,
        sort_by=arxiv.SortCriterion.Relevance,
        sort_order=arxiv.SortOrder.Descending,
    )

    # Download 
    downloaded = 0
    seen = 0
    while True:
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
                    print(f"[{downloaded}/{NUM_PAPERS}] saved {arxiv_id}")
                except Exception as e:
                    print(f"[skip] {arxiv_id}: {e}", file=sys.stderr)
            break
        except arxiv.UnexpectedEmptyPageError as e:
            # Rare hiccup; quick retry
            print("[warn] empty page from API; retrying in 5s…", file=sys.stderr)
            time.sleep(5)
            continue

print(f"[done] seen={seen}, newly downloaded={downloaded}, saved to {OUT_DIR}")
