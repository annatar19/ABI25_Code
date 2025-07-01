#!/bin/bash

####### User‐configurable macros #######
# Number of times to run each benchmark:
RUNS=150

# Cache mode: 0 = cold cache, 1 = hot cache, 2 = both
CACHE_FLAG=2

# Output mode: 0 = scatter plot, 1 = mean, std and throughput
STAT_FLAG=0

# Directories (relative to this script)
INPUT_DIR="./input"
BUILD_DIR="./build"
BENCH="$BUILD_DIR/benchmarks"
########################################set -e

mkdir -p "$INPUT_DIR"

rm *.csv 2>/dev/null

# This creates soft links for the scaled outputs and the converted outputs. It does assume that those exists, if not you need to run those run.sh's first.
cd "$INPUT_DIR" || exit
# The RNTuple files
source_dir=../../scaling/output
for src in "$source_dir"/*B2HHH*; do
    dest="$(basename "$src")"
    ln -s "$src" "$dest" 2>/dev/null || :
done
# The ORC files
source_dir=../../rntuple_to_orc/output
for src in "$source_dir"/*B2HHH*; do
    dest="$(basename "$src")"
    ln -s "$src" "$dest" 2>/dev/null || :
done
cd ../

if [[ ! -x "$BENCH" ]]; then
    echo "Error: benchmark binary not found or not executable at $BENCH" >&2
    exit 1
fi

shopt -s nullglob
# todo if scatter then 1x
for file in "$INPUT_DIR"/1x_*; do
    case "$file" in
    *.ntuple.root)
        FORMAT_FLAG=0
        ;;
    *.orc)
        FORMAT_FLAG=1
        ;;
    *)
        echo "Skipping unknown file type: $file"
        continue
        ;;
    esac

    if [[ "$CACHE_FLAG" -eq 2 ]]; then
        CACHE_MODES=(0 1)
    else
        CACHE_MODES=($CACHE_FLAG)
    fi

    # echo "Running: $BENCH $file $RUNS $FORMAT_FLAG $CACHE_FLAG $STAT_FLAG"
    # "$BENCH" "$file" "$RUNS" "$FORMAT_FLAG" "$CACHE_FLAG" "$STAT_FLAG"
    for mode in "${CACHE_MODES[@]}"; do
        cache_desc=$([[ "$mode" -eq 0 ]] && echo "cold" || echo "hot")
        echo "Running ($cache_desc cache): $BENCH $file $RUNS $FORMAT_FLAG $mode $STAT_FLAG"
        "$BENCH" "$file" "$RUNS" "$FORMAT_FLAG" "$mode" "$STAT_FLAG"
    done
done
