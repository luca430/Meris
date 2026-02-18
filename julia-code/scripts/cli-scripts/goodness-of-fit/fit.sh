#!/bin/bash
set -e

# Get the directory of the script, irregardless of where it is called from
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# julia --project gutenberg.jl
# julia --project arxiv.jl

julia --project gtex.jl
