# simPed
a pedigree simulation for complex and Mendelian traits

## Setup
This is a list of required R packages: 
- tidyverse: manipulation of dataframes
- MASS: to sample from multvariate normal distribution
- mvnfast: to sample from multvariate normal distribution, fast
- pedtools: generates pedigree structures, for plotting and calculating kinship matrix
- ribd: calculate kinship matrix for pedigrees
- RColorBrewer: for plotting polygenic G values
```
packages = c("tidyverse", "MASS", "mvnfast", "pedtools", "ribd", "RcolorBrewer")
install.packages(packages)
```

