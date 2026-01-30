#!/usr/bin/env python3
import argparse
import os
import sys
import tarfile
import zipfile
from urllib.request import Request, urlopen

from matplotlib.pylab import extract

URLS = {
    ("txt", "tar"): "https://www.rfc-editor.org/in-notes/tar/RFC-all.tar.gz",
    ("txt", "zip"): "https://www.rfc-editor.org/in-notes/tar/RFC-all.zip",
    # You can also use the bulk page to grab PDF/XML archives. :contentReference[oaicite:3]{index=3}
}

def download(url: str, out_path: str, chunk_size: int = 1024 * 1024) -> None:
    req = Request(url, headers={"User-Agent": "rfc-bulk-downloader/1.0"})
    with urlopen(req) as r, open(out_path, "wb") as f:
        total = r.headers.get("Content-Length")
        total = int(total) if total else None
        done = 0
        while True:
            chunk = r.read(chunk_size)
            if not chunk:
                break
            f.write(chunk)
            done += len(chunk)
            if total:
                pct = done * 100.0 / total
                print(f"\rDownloaded: {done}/{total} bytes ({pct:5.1f}%)", end="", file=sys.stderr)
            else:
                print(f"\rDownloaded: {done} bytes", end="", file=sys.stderr)
    print(file=sys.stderr)

def extract_txt_only(archive_path: str, out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)

    if archive_path.endswith((".tar.gz", ".tgz")):
        import tarfile
        with tarfile.open(archive_path, "r:gz") as t:
            members = [m for m in t.getmembers()
                       if m.isfile() and m.name.lower().endswith(".txt")]
            t.extractall(path=out_dir, members=members)

    elif archive_path.endswith(".zip"):
        import zipfile
        with zipfile.ZipFile(archive_path) as z:
            for name in z.namelist():
                if name.lower().endswith(".txt") and not name.endswith("/"):
                    z.extract(name, path=out_dir)
    else:
        raise ValueError(f"Unknown archive type: {archive_path}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fmt", choices=["txt"], default="txt", help="Archive content format (default: txt)")
    ap.add_argument("--kind", choices=["tar", "zip"], default="tar", help="Archive type (default: tar)")
    ap.add_argument("--out", default="raw-data", help="Output directory (default: raw-data)")
    ap.add_argument("--keep-archive", action="store_true", help="Keep downloaded archive file")
    args = ap.parse_args()

    url = URLS.get((args.fmt, args.kind))
    if not url:
        raise SystemExit("No URL configured for that (fmt, kind). Use --fmt txt --kind tar|zip.")

    os.makedirs(args.out, exist_ok=True)
    archive_name = os.path.basename(url)
    archive_path = os.path.join(args.out, archive_name)

    print(f"Downloading {url}", file=sys.stderr)
    download(url, archive_path)

    print(f"Extracting to {args.out}", file=sys.stderr)
    extract_txt_only(archive_path, args.out)

    if not args.keep_archive:
        os.remove(archive_path)

    print("Done.", file=sys.stderr)

if __name__ == "__main__":
    main()
