#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SRC_DIR="./raw-data"
OUT_DIR="./processed"
TMP_DIR="./tmp"

mkdir -p "$OUT_DIR" "$TMP_DIR"
shopt -s nullglob

# Recursively find all archives
find "$SRC_DIR" -type f \( -name "*.tar.gz" -o -name "*.tgz" \) | sort | while IFS= read -r f; do
  # Continue on errors instead of aborting
  {
    rel_path="${f#$SRC_DIR/}"
    base_name="$(basename "$f" .tar.gz)"
    base_name="${base_name%.tgz}"

    out_subdir="$(dirname "$rel_path")"
    out_dir="$OUT_DIR/$out_subdir"
    mkdir -p "$out_dir"
    final="$out_dir/$base_name.txt"

    [[ -s "$final" ]] && { echo "✔ $final"; continue; }

    rm -rf "$TMP_DIR"/*; mkdir -p "$TMP_DIR"
    if ! tar -xzf "$f" -C "$TMP_DIR" 2>/dev/null; then
      echo "[warn] cannot extract $f"
      continue
    fi

    if ! find "$TMP_DIR" -type f -name "*.tex" | grep -q .; then
      : > "$final"
      echo "[warn] no .tex in $f"
      continue
    fi

    find "$TMP_DIR" -type f -name "*.tex" | sort | xargs cat \
      | awk '/\\begin{document}/, /\\end{document}/' \
      | detex -w -n \
      | tr '[:upper:]' '[:lower:]' \
      | perl -0777 -pe '
          s/\r//g;
          s/([[:alpha:]])-\s*\n\s*([[:alpha:]])/$1$2/g;
          s/^[[:space:]]*-\s+([[:alpha:]])/$1/gm;
          s/([[:space:]])-\s*([[:alpha:]])/$1$2/g;
          s/[[:space:]]?[–—][[:space:]]?/ /g;
          s/[ \t]{2,}/ /g;
        ' \
      | grep -viE '^[[:space:]]*([0-9ivx\.]*[[:space:]]*)?(introduction|results?|discussion|methods?|simulations?|conclusion|abstract|acknowledg(e)?ments?)\b' \
      > "$final"

    echo "✔ wrote $final"
  } || echo "[error] failed on $f"
done

rm -rf "$TMP_DIR"
echo "✅ All done. Clean text files are in: $OUT_DIR"
