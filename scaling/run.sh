#!/bin/bash

set -e

FILES_LIST="files.txt"
SCALE_UPPER_LIMIT=1

mkdir -p output

while IFS=' ' read -r first second; do
    cd build/
    for i in 0.0001 0.001 0.01 0.1 0.25 0.5 0.75; do
        echo "Scaling $first to $i"
        ./rntuple_scaling "$first" "$second" "$i" >/dev/null
    done
    for ((i = 1; i <= SCALE_UPPER_LIMIT; i++)); do
        echo "Scaling $first to $i"
        ./rntuple_scaling "$first" "$second" "$i" >/dev/null
    done
    cd ../
done <"$FILES_LIST"
