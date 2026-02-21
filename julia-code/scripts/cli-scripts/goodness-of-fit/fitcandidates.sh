#!/bin/bash
set -e

# Get the directory of the script, irregardless of where it is called from
SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
NEPSILON=100

for file in "$SCRIPTDIR"/*.jl; do
    [ -f "$file" ] || continue
    julia --project "$file" -n $NEPSILON
done


