#!/bin/bash
set -e

# Get the directory of the script, irregardless of where it is called from
SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
NEPSILON=100
NTOP=50
# tip: check `julia --project [name].jl --help` for arg options
julia --threads 8 --project arxiv.jl -f -n $NEPSILON --top $NTOP
# julia --project finance.jl -f -n $NEPSILON -t $NTOP
# julia --project gowalla.jl -f -n $NEPSILON -t $NTOP
julia --project gutenberg.jl -f -n $NEPSILON -t $NTOP
julia --project lego.jl -f -n $NEPSILON -t $NTOP
julia --project rfc.jl -f -n $NEPSILON -t $NTOP
julia --project bcitree.jl -f -n $NEPSILON -t $NTOP
julia --project biotime.jl -f -n $NEPSILON -t $NTOP
julia --project gtex.jl -f -n $NEPSILON -t $NTOP
julia --project otu.jl -f -n $NEPSILON -t $NTOP


