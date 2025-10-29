#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

SRC_DIR="./Data/arx_sources/nlin/AO/sources"
TXT_DIR="./Data/arx_sources/nlin/AO/texts"
TMP_DIR="./tmp"

mkdir -p "$TXT_DIR" "$TMP_DIR"
shopt -s nullglob

for f in "$SRC_DIR"/*.tar.gz "$SRC_DIR"/*.tgz; do
  name="$(basename "$f" .tar.gz)"; name="${name%.tgz}"
  final="$TXT_DIR/$name.txt"

  [[ -s "$final" ]] && { echo "✔ $final"; continue; }

  rm -rf "$TMP_DIR"/*; mkdir -p "$TMP_DIR"
  if ! tar -xzf "$f" -C "$TMP_DIR"; then
    echo "[warn] cannot extract $f"; continue
  fi

  if ! find "$TMP_DIR" -type f -name "*.tex" | grep -q .; then
    : > "$final"; echo "[warn] no .tex in $f"; continue
  fi

  # Strip LaTeX → lowercase → de-hyphenate & clean bullets/dashes → drop section headers
  find "$TMP_DIR" -type f -name "*.tex" | sort | xargs cat \
    | awk '/\\begin{document}/, /\\end{document}/' \
    | detex -w -n \
    | tr '[:upper:]' '[:lower:]' \
    | perl -0777 -pe '
        s/\r//g;

        # 1) Join words split by hyphenation at line breaks: "galax-\nies" → "galaxies"
        s/([[:alpha:]])-\s*\n\s*([[:alpha:]])/$1$2/g;

        # 2) Remove list/bullet hyphens at start of lines: "- galaxies" → "galaxies"
        s/^[[:space:]]*-\s+([[:alpha:]])/$1/gm;

        # 3) Remove stray hyphen directly before a word within a line: " -galaxies" → " galaxies"
        s/([[:space:]])-\s*([[:alpha:]])/$1$2/g;

        # (Optional) normalize en/em dashes to spaces
        s/[[:space:]]?[–—][[:space:]]?/ /g;

        # 4) Squeeze repeated spaces
        s/[ \t]{2,}/ /g;
      ' \
    | grep -viE '^[[:space:]]*([0-9ivx\.]*[[:space:]]*)?(introduction|results?|discussion|methods?|simulations?|conclusion|abstract|acknowledg(e)?ments?)\b' \
    > "$final"

  echo "✔ wrote $final"
done

rm -rf "$TMP_DIR"
echo "✅ All done. Clean text files are in: $TXT_DIR"
