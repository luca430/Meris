#!/usr/bin/env bash
set -euo pipefail

RFCURL="https://www.rfc-editor.org/in-notes/tar/RFC-all.tar.gz"
TARCHIVE="RFC-all.tar.gz"
OUT="raw-data"

# Download and unpack
if [[ -f "$TARCHIVE" ]]; then
    echo "Tarball already exists: $TARCHIVE"
else
    curl -fL "$RFCURL" -o "$TARCHIVE"
fi
mkdir -p "$OUT"
tar --verbose --extract --file="$TARCHIVE" --directory="$OUT" --wildcards 'rfc*.txt'
# Parse
# [...]
