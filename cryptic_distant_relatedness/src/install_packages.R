required_packages <- c(
  "argparse", "data.table", "R.utils", "tidyverse", "MASS", "mvnfast", "kinship2", "sn", "igraph",
  "pedtools", "ribd", "RColorBrewer", "this.path", "wesanderson", "broom", "fmsb", "patchwork", "dplyr", "tidyr"
)

for (package in required_packages) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, repos = "https://cloud.r-project.org")
  }
}

# Create an empty file to indicate successful installation
file.create("logs/packages_installed.txt")