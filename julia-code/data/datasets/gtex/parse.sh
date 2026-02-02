#!/usr/bin/env bash
set -euo pipefail

HEADER_ROW=3
RESUME=0

usage() {
  cat <<EOF
Usage: $0 -i RAW_ROOT -o OUT_ROOT [-r HEADER_ROW] [--resume]

Input layout:  RAW_ROOT/<tissue>/*.gz
Output layout: OUT_ROOT/<tissue>/*.nonzero_cols.gz + *.nonzero_samples.txt

  -i  Root directory (e.g. raw-data)
  -o  Output root directory (e.g. processed)
  -r  Header row (1-based). Default: 3
  --resume  Skip file if output already exists

EOF
}

RAW_ROOT=""
OUT_ROOT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i) RAW_ROOT="$2"; shift 2 ;;
    -o) OUT_ROOT="$2"; shift 2 ;;
    -r) HEADER_ROW="$2"; shift 2 ;;
    --resume) RESUME=1; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "$RAW_ROOT" || -z "$OUT_ROOT" ]]; then
  echo "Error: -i and -o are required" >&2
  usage
  exit 2
fi
if [[ ! -d "$RAW_ROOT" ]]; then
  echo "Error: input directory not found: $RAW_ROOT" >&2
  exit 2
fi
mkdir -p "$OUT_ROOT"

# zcat or gzcat
ZCAT_CMD="zcat"
if ! command -v zcat >/dev/null 2>&1; then
  if command -v gzcat >/dev/null 2>&1; then
    ZCAT_CMD="gzcat"
  else
    echo "Error: need zcat or gzcat in PATH" >&2
    exit 2
  fi
fi

# REQUIRE gawk (mawk often truncates wide GTEx rows)
if ! command -v gawk >/dev/null 2>&1; then
  echo "Error: gawk not found. Install it: sudo apt-get install -y gawk" >&2
  exit 2
fi
AWK="gawk"

if ! command -v gzip >/dev/null 2>&1; then
  echo "Error: need gzip in PATH" >&2
  exit 2
fi

mapfile -t FILES < <(find "$RAW_ROOT" -mindepth 2 -maxdepth 2 -type f -name "*.gz" | sort)
if [[ "${#FILES[@]}" -eq 0 ]]; then
  echo "No .gz files found under: $RAW_ROOT/<tissue>/*.gz" >&2
  exit 0
fi

echo "[info] using: $AWK"
echo "[info] files: ${#FILES[@]}"

for INFILE in "${FILES[@]}"; do
  TISSUE="$(basename "$(dirname "$INFILE")")"
  BASENAME="$(basename "$INFILE")"

  OUT_DIR="$OUT_ROOT/$TISSUE"
  mkdir -p "$OUT_DIR"

  OUTFILE="$OUT_DIR/${BASENAME%.gz}.nonzero_cols.gz"
  KEEPFILE="$OUT_DIR/${BASENAME%.gz}.nonzero_samples.txt"

  if [[ "$RESUME" -eq 1 && -f "$OUTFILE" && -s "$OUTFILE" ]]; then
    echo "[skip] $TISSUE / $BASENAME"
    continue
  fi

  echo "[*] $TISSUE / $BASENAME"

  # Pass 1: detect nonzero columns (iterate to header width maxNF)
  "$ZCAT_CMD" "$INFILE" | \
  "$AWK" -F'\t' -v HR="$HEADER_ROW" '
    NR==HR { for(i=2;i<=NF;i++) name[i]=$i; maxNF=NF; next }
    NR>HR  {
      for(i=2;i<=maxNF;i++){
        v = (i<=NF ? $i : "")
        if(!nz[i] && v!="" && v!="0") nz[i]=1
      }
    }
    END { for(i=2;i<=maxNF;i++) if(nz[i]) print name[i] }
  ' > "$KEEPFILE"

  KEPT_N="$(wc -l < "$KEEPFILE" | tr -d " ")"
  echo "  [ok] nonzero columns: $KEPT_N"

  # Pass 2: filter file, padding missing fields up to maxNF
  "$ZCAT_CMD" "$INFILE" | \
  "$AWK" -F'\t' -v HR="$HEADER_ROW" -v KEEPFILE="$KEEPFILE" '
    BEGIN {
      while ((getline line < KEEPFILE) > 0) keepname[line]=1
      close(KEEPFILE)
    }
    NR < HR { print; next }
    NR == HR {
      maxNF = NF
      for(i=1;i<=maxNF;i++) if(i==1 || keepname[$i]) keepi[i]=1

      first=1
      for(i=1;i<=maxNF;i++) if(keepi[i]) {
        printf "%s%s", (first?"":"\t"), $i
        first=0
      }
      printf "\n"
      next
    }
    NR > HR {
      first=1
      for(i=1;i<=maxNF;i++) if(keepi[i]) {
        v = (i<=NF ? $i : "")
        printf "%s%s", (first?"":"\t"), v
        first=0
      }
      printf "\n"
    }
  ' | gzip -c > "$OUTFILE"

  echo "  [saved] $OUTFILE"
done

echo "[done]"
