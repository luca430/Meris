#!/usr/bin/env bash

SRC_DIR="./Data/arx_sources/astro-ph/CO/sources"
TXT_DIR="./Data/arx_sources/astro-ph/CO/texts"
TMP_DIR="./tmp"

mkdir -p "$TXT_DIR" "$TMP_DIR"
shopt -s nullglob

for f in "$SRC_DIR"/*.tar.gz "$SRC_DIR"/*.tgz; do
  name="$(basename "$f" .tar.gz)"; name="${name%.tgz}"
  final="$TXT_DIR/$name.txt"

  # Skip if already done
  [[ -s "$final" ]] && { echo "✔ $final"; continue; }

  rm -rf "$TMP_DIR"/*; mkdir -p "$TMP_DIR"
  if ! tar -xzf "$f" -C "$TMP_DIR"; then
    echo "[warn] cannot extract $f"; continue
  fi

  if ! find "$TMP_DIR" -type f -name "*.tex" | grep -q .; then
    : > "$final"; echo "[warn] no .tex in $f"; continue
  fi

  # Convert all .tex -> plain, then scrub math, then normalize to words
  find "$TMP_DIR" -type f -name "*.tex" -exec pandoc -f latex -t plain --strip-comments {} + \
  | perl -0777 -pe '
      s/\$\$(.|\n)*?\$\$//g;
      s/\\\[(.|\n)*?\\\]//g;
      s/\$(.|\n)*?\$//g;
      s/\\\((.|\n)*?\\\)//g;
      for my $env (qw(equation align gather multline eqnarray flalign aligned split cases)) {
        s/\\begin\{$env\*?\}(.|\n)*?\\end\{$env\*?\}//g;
      }
      s/\\[a-zA-Z@]+(\s*\[[^\]]*\])?(\s*\{[^}]*\})?//g;
      s/[{}]//g;
    ' \
  | awk 'BEGIN{IGNORECASE=1} /^references[[:space:]]*$/ {exit} {print}' \
  | tr '[:upper:]' '[:lower:]' | sed 's/[^a-z]/ /g' | tr -s " " \
  > "$final"

  echo "✔ wrote $final"
done
