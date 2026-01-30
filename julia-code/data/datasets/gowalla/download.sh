#!/usr/bin/env bash
set -euo pipefail

DIR="${1:-raw-data}"
mkdir -p "$DIR"

BASE="https://snap.stanford.edu/data"

FILES=(
  "loc-gowalla_edges.txt.gz"
  "loc-gowalla_totalCheckins.txt.gz"
)

echo "Downloading Gowalla dataset to $DIR"
for f in "${FILES[@]}"; do
  echo "Downloading $f..."
  wget -c "$BASE/$f" -O "$DIR/$f"
done

echo "Done! Files in $DIR:"
ls -lh "$DIR"
