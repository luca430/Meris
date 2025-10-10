# Script to download the newest arXiv source tarballs for a given category.
import os, time, sys
import arxiv

# --- config ---
CATEGORY      = "nlin.SI"           # e.g. "astro-ph.CO"
NUM_PAPERS    = 200
DELAY_SECONDS = 3                        # be polite (≥3)
OUT_DIR       = os.path.join(os.path.dirname(__file__),
                             "Data", f"arx_sources/{CATEGORY[:-3]}/{CATEGORY[-2:]}/sources")
os.makedirs(OUT_DIR, exist_ok=True)

search = arxiv.Search(
    query=f"cat:{CATEGORY}",
    max_results=NUM_PAPERS,              # still set this, but we'll hard-cap ourselves too
    sort_by=arxiv.SortCriterion.SubmittedDate,
    sort_order=arxiv.SortOrder.Descending,
)

client = arxiv.Client(delay_seconds=DELAY_SECONDS)

# Seed seen IDs with what's already on disk (so we don't redownload or recount them)
ids_seen = {
    fn[:-7]  # strip ".tar.gz"
    for fn in os.listdir(OUT_DIR)
    if fn.endswith(".tar.gz")
}

downloaded = 0        # newly downloaded this run
processed_unique = 0  # unique results considered this run (including already-on-disk)

max_retries = 4         # stop after this many consecutive empty pages
retry_count = 0
backoff_base = 5        # seconds; grows exponentially

while True:
    try:
        made_progress = False
        for r in client.results(search):
            arx_id = r.entry_id.rsplit("/", 1)[-1]

            # de-dup & hard cap (assume ids_seen and processed_unique exist)
            if arx_id in ids_seen:
                continue
            ids_seen.add(arx_id)
            processed_unique += 1
            made_progress = True

            if processed_unique > NUM_PAPERS:
                break

            out_path = os.path.join(OUT_DIR, f"{arx_id}.tar.gz")
            if not os.path.exists(out_path):
                arx_id = r.entry_id.rsplit("/", 1)[-1]
                short_name = f"{arx_id}.tar.gz"
                r.download_source(dirpath=OUT_DIR, filename=short_name)
                downloaded += 1
                print(f"[{min(processed_unique, NUM_PAPERS)}/{NUM_PAPERS}] saved {arx_id}")

        # finished the iterator or hit the cap → done
        break

    except arxiv.UnexpectedEmptyPageError:
        # If we already hit the cap, stop anyway
        if processed_unique >= NUM_PAPERS:
            break

        retry_count += 1
        if retry_count > max_retries:
            print(f"[warn] too many empty pages; stopping after {max_retries} retries.", file=sys.stderr)
            break

        # Exponential backoff (5, 10, 20, …)
        sleep_s = backoff_base * (2 ** (retry_count - 1))
        print(f"[warn] empty page from API; retrying in {sleep_s}s… (attempt {retry_count}/{max_retries})", file=sys.stderr)
        time.sleep(sleep_s)
        continue


print(f"[done] processed_unique={processed_unique}, newly_downloaded={downloaded}, saved to {OUT_DIR}")
