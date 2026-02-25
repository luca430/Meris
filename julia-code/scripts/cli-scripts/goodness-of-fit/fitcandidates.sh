#!/bin/bash
set -e

# Get the directory of the script, irregardless of where it is called from
SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
NEPSILON=64
NTOP=10
# tip: check `julia --project [name].jl --help` for arg options
julia --project arxiv.jl -f -n $NEPSILON -t $NTOP
julia --project finance.jl -f -n $NEPSILON -t $NTOP
julia --project gowalla.jl -f -n $NEPSILON -t $NTOP
julia --project gutenberg.jl -f -n $NEPSILON -t $NTOP
julia --project lego.jl -f -n $NEPSILON -t $NTOP
julia --project rfc.jl -f -n $NEPSILON -t $NTOP


