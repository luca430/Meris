#!/usr/bin/env bash
shopt -s nullglob

DIRECTORY="raw-text/"
OUTCSV="rfc-documentsize.csv"

echo "sample_id,documentsize" > $OUTCSV

for file in "$DIRECTORY"/rfc[0-9]*.txt; do
    #~ Extract the name using `basename` and `cut`
    name=$(basename "$file" | cut -d. -f1)
    #~ Count words
    count=$(wc -w < "$file")
    #~ echo name and counts to OUTCSV
    echo "$name,$count" >> $OUTCSV
done
