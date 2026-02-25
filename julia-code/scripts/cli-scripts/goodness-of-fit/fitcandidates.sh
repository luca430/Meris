#!/bin/bash
set -e

# Get the directory of the script, irregardless of where it is called from
SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
NEPSILON=100

# julia --project otu.jl -n $NEPSILON
# julia --project rfc.jl -n $NEPSILON

for file in "$SCRIPTDIR"/*.jl; do
    [ -f "$file" ] || continue
    echo "$file"
    julia --project "$file" -n $NEPSILON
done


