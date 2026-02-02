#!/usr/bin/env python3
import argparse
import os
import re
import time
from pathlib import Path
from datetime import datetime
import requests
import json

GUTENDEX = "https://gutendex.com/books"
UA = "luca-gutendex-downloader/1.1"

def safe_name(s: str) -> str:
    s = re.sub(r"[^\w\s\-.(),]", "", s, flags=re.UNICODE).strip()
    s = re.sub(r"\s+", "_", s)
    return s[:140] if s else "book"

def pick_format(formats: dict, prefer: str):
    if prefer == "txt":
        candidates = [k for k in formats.keys() if k.startswith("text/plain")]
        candidates += [k for k in formats.keys() if k.startswith("text/html")]
    else:
        candidates = [k for k in formats.keys() if k.startswith("application/epub+zip")]

    for mime in candidates:
        url = formats.get(mime)
        if isinstance(url, str) and url.startswith("http"):
            return mime, url
    return None

def ext_from_mime(mime: str) -> str:
    if mime.startswith("application/epub+zip"):
        return ".epub"
    if mime.startswith("text/plain"):
        return ".txt"
    if mime.startswith("text/html"):
        return ".html"
    return ".bin"

def iter_gutendex(params: dict, session: requests.Session, sleep_s: float):
    url = GUTENDEX
    first = True
    while url:
        r = session.get(url, params=params if first else None, timeout=60)
        r.raise_for_status()
        data = r.json()
        yield from data.get("results", [])
        url = data.get("next")
        first = False
        time.sleep(sleep_s)

def download_file(url: str, out_path: Path, session: requests.Session, sleep_s: float) -> bool:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if out_path.exists():
        return False

    with session.get(url, stream=True, timeout=120) as r:
        r.raise_for_status()
        tmp = out_path.with_suffix(out_path.suffix + ".part")
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 256):
                if chunk:
                    f.write(chunk)
        tmp.replace(out_path)

    time.sleep(sleep_s)
    return True

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=50, help="Max books to download per language.")
    ap.add_argument("--langs", default="en", help="Comma-separated language codes, e.g. en,it,fr")
    ap.add_argument("--q", default="", help="Search query (optional).")
    ap.add_argument("--prefer", choices=["txt", "epub"], default="txt", help="Preferred format.")
    ap.add_argument("--sort", choices=["ascending", "descending", "popular"], default="descending",
                    help="Sorting on Gutendex (use ascending/descending for replicability).")
    ap.add_argument("--out", default="raw-txt", help="Output directory.")
    ap.add_argument("--sleep", type=float, default=1.0, help="Seconds between requests.")
    ap.add_argument("--manifest", default="", help="Optional manifest file path (jsonl).")
    args = ap.parse_args()

    # You said you want replicable: popular is NOT replicable.
    if args.sort == "popular":
        print("You are wrong to use --sort popular if you want replicability. Use ascending or descending.", flush=True)

    langs = [x.strip() for x in args.langs.split(",") if x.strip()]
    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({"User-Agent": UA})

    manifest_path = Path(args.manifest) if args.manifest else None
    mf = open(manifest_path, "a", encoding="utf-8") if manifest_path else None

    try:
        for lang in langs:
            lang_dir = out_root / lang
            lang_dir.mkdir(parents=True, exist_ok=True)

            params = {"languages": lang, "sort": args.sort}
            if args.q:
                params["search"] = args.q

            downloaded = 0
            seen = 0

            for book in iter_gutendex(params, session, args.sleep):
                if downloaded >= args.n:
                    break

                seen += 1
                book_id = book.get("id")
                title = book.get("title", f"book_{book_id}")
                formats = book.get("formats", {})

                picked = pick_format(formats, args.prefer)
                if not picked:
                    continue
                mime, url = picked

                fname = f"{book_id:06d}_{safe_name(title)}{ext_from_mime(mime)}"
                out_path = lang_dir / fname

                try:
                    did = download_file(url, out_path, session, args.sleep)
                    if did:
                        downloaded += 1
                        print(f"[{lang}] [{downloaded}/{args.n}] {title} -> {out_path}")
                
                        if mf:
                            record = {
                                "ts": datetime.utcnow().isoformat() + "Z",
                                "lang": lang,
                                "sort": args.sort,
                                "q": args.q,
                                "id": book_id,
                                "title": title,
                                "mime": mime,
                                "url": url,
                            }
                            mf.write(json.dumps(record, ensure_ascii=False) + "\n")
                            mf.flush()
                
                except Exception as e:
                    print(f"[skip] lang={lang} id={book_id} {title}: {e}")

            print(f"[done] lang={lang} seen={seen}, downloaded={downloaded}, saved to {lang_dir}")
    finally:
        if mf:
            mf.close()

if __name__ == "__main__":
    main()
