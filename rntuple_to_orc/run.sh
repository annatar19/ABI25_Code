#!/bin/bash

set -e
mkdir -p output input

cd input
source_dir=../../scaling/output
for src in "$source_dir"/*; do
    dest="$(basename "$src")"
    ln -s "$src" "$dest" 2>/dev/null || :
done
cd ../

while IFS=' ' read -r postfix access; do
    for file in input/*"$postfix"; do
        echo "Converting $file to ORC…"

        cd build
        ./rntuple_to_orc "../$file" "$access" >/dev/null
        cd ..
    done
done <files.txt
