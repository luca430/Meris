#!/usr/bin/env bash
set -euo pipefail

DIR="${1:-raw-data}"
mkdir -p "$DIR"

BASE="https://snap.stanford.edu/data"

FILES=(
  "loc-brightkite_edges.txt.gz"
  "loc-brightkite_totalCheckins.txt.gz"
)

echo "Downloading Brightkite dataset to $DIR"
for f in "${FILES[@]}"; do
  echo "Downloading $f..."
  wget -c "$BASE/$f" -O "$DIR/$f"
done

echo "Done! Files in $DIR:"
ls -lh "$DIR"