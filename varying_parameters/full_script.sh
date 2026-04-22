#!/bin/bash
for dir in */; do
    if [[ "${dir%/}" == "src" ]]; then continue; fi
    echo "Running snakemake in ${dir%/}..."
    (cd "$dir" && \
     ls -al results/plots/* && \
     cd ..)
done