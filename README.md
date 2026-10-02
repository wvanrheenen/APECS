# APECS
A framework to simulate ALS and ALS-associated disease under a monogenic/Mendelian and polygenic/complex disease model.

## Github Architecture
In this main directory, you'll find 'snakefile', the pipeline needed to run the simulations, as per our main analyses in the article. The results of this analysis will be written to the 'results' subdirectory. 

This is a list of the subdirectories in our main github page:
- interactive_tutorial:
  - Here, you'll find an interactive R script with an introductory explanation file
- src:
  - This directory contains the main functions required to run all analyses. 
  - Further subdirectories contain adaptations of the main functions found in this src directory. 
- varying_parameters:
  - This directory contains the files needed for the sensitivity analyses and demographic changes analyses in our main article.
- rationale:
  - This directory contains the data required to validate whether the simulations match input parameters, real world demographic data and if the simulations actually simulate phenotypic traits as is expected.
  - A separate README.md file is available for explanation on each rationale step.
- data:
  - This directory contains gathered and processed real world demographic and genetic data required as input to run the simulations.
- ALPINE_plug:
  - This directory contains a pipeline needed to run analyses on cryptic distant relatives in monogenic index patient pedigrees.
- rshiny_slider:
  - This directory contains the pipeline used to generate the 'Variable parameters' tab of the web-based probability calculator

## Setup
### Snakefile
**Snakemake + SLURM executor plugin** required for HPC execution.

#### Installation (Conda - Recommended)
```bash
# Install snakemake + SLURM executor plugin
conda create -n APECS_snakemake -c conda-forge -c bioconda snakemake snakemake-executor-plugin-slurm
conda activate APECS_snakemake
```

#### Run Pipeline (Login Node Only)
```bash
conda activate APECS_snakemake

# This commnand can be applied for all snakefiles in this github; change the -s {snakefile_name} to the correct name where necessary for parallelization
snakemake -s snakefile --executor slurm --jobs 500 \
  --default-resources mem_mb=1000 runtime=900 constraint="avx2" \
  --rerun-incomplete --keep-going --latency-wait 90 --cores 1 \
  --slurm-keep-successful-logs all
```


### R Packages
This is a list of required R packages:
- argparse: command line argument parsing
- data.table: fast data manipulation
- R.utils: utility functions
- tidyverse: manipulation of dataframes (dplyr, tidyr, etc.)
- MASS: to sample from multivariate normal distribution
- mvnfast: to sample from multivariate normal distribution, fast
- kinship2: pedigree and kinship analysis
- sn: skew-normal distribution sampling
- igraph: network analysis and visualization
- pedtools: generates pedigree structures, for plotting and calculating kinship matrix
- ribd: calculate kinship matrix for pedigrees
- RColorBrewer: for plotting polygenic G values
- this.path: get script path
- wesanderson: color palettes for plotting in article
- broom: convert statistical objects to tidy data
- fmsb: miscellaneous functions (radar charts)
- patchwork: plot composition
- dplyr: data manipulation (tidyverse core)
- tidyr: data tidying (tidyverse core)

```r
required_packages <- c(
  "argparse", "data.table", "R.utils", "tidyverse", "MASS", "mvnfast", 
  "kinship2", "sn", "igraph", "pedtools", "ribd", "RColorBrewer", 
  "this.path", "wesanderson", "broom", "fmsb", "patchwork", "dplyr", "tidyr"
)
install.packages(required_packages)
```

These are loaded via `src/libraries_simPed.R`.