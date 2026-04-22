#!/bin/bash

# Step 1: create simulations
snakemake -s snakefile --executor slurm --jobs 500 --default-resources mem_mb=5000 constraint="avx2" --rerun-incomplete --keep-going --latency-wait 90 --cores 4 --slurm-keep-successful-logs --forceall all

# Step 2: concatenate simulations
cd results

for base_dir in h2_als K_als; do
  echo "Processing $base_dir..."

  # Loop over each subdirectory inside the base_dir
  for subdir in "$base_dir"/*; do
    if [ -d "$subdir" ]; then
      param=$(basename "$subdir")  # e.g. varying_gen
      output_file="$base_dir/${param}.all"

      echo "Concatenating all .out files in $subdir into $output_file"

      shopt -s nullglob
      files=("$subdir"/*.out)

      if [ ${#files[@]} -gt 0 ]; then
        head -n 1 "${files[0]}" > "$output_file"
        tail -n +2 -q "${files[@]}" >> "$output_file"
      else
        echo "No .out files found in $subdir, skipping."
      fi
    fi
  done
done

cd ..

# Step 3: create visualization
Rscript polygenic_visualization.R