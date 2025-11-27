#!/usr/bin/bash

# Specify no. of papers [default=10]
NUM_PAPERS=${1:-100}
# Fetch NUM_PAPERS papers from the categories defined in `categories.txt`
# note: assumes `poetry`
# note: if not using `poetry`, change here (e.g., `uv`, or regular `python`)
poetry run python query-arxiv.py --n $NUM_PAPERS

# Prevent unmatched globs (e.g. *.txt) from expanding to themselves;
# make them expand to nothing instead
shopt -s nullglob

# Define temporary dir. for unzip/untar
TMP_DIR="./tmp"

while IFS= read -r line; do
    IFS="." read -r CAT SUBCAT <<< "$line"
    SRC_DIR="./raw-data/$CAT/$SUBCAT"
    TXT_DIR="./raw-text/$CAT/$SUBCAT"
    mkdir -p $TXT_DIR
    echo "Checking... $SRC_DIR"

    for f in "$SRC_DIR"/*.tar.gz "$SRC_DIR"/*.tgz; do
        [ -e "$f" ] || continue   # skip if glob did not match
        # Get name
        name="$(basename "$f" .tar.gz)"; name="${name%.tgz}"
        final="$TXT_DIR/$name.txt"
        # echo "Name: $name"

        # Extract
        rm -rf "$TMP_DIR"/*; mkdir -p "$TMP_DIR"
        if ! tar -xzf "$f" -C "$TMP_DIR"; then
            # echo "[warn] cannot extract $f"; continue
            continue
        fi
        
        if ! find "$TMP_DIR" -type f -name "*.tex" | grep -q .; then
            # : > "$final"; echo "[warn] no .tex in $f"; continue
            continue
        fi
        
        # Use `detex` to strip all LaTeX
        # Also convert to lower case using `tr`, and finally
        # `grep` to get only lines with alphabetical characters
        find "$TMP_DIR" -type f -name "*.tex" -print0 | sort -z | xargs -0 cat \
             |awk '
                  /\\begin{document}/ {inrange=1; next}
                  /\\section\*?\{[Aa]cknowledg/ {inrange=0}
                  /\\begin{thebibliography}/ {inrange=0}
                  /\\end{document}/ {inrange=0}
                  inrange
                  ' \
            | detex -w -n \
            | tr '[:upper:]' '[:lower:]' \
            | grep -E '^[a-z]+$' > "$final"
        # If detex produced nothing, remove or warn
        if [ ! -s "$final" ]; then
            echo "[warn] empty output for $name — removing"
            rm -f "$final"
        else
            echo "✔ wrote $final"
        fi
    done
done < categories.txt
