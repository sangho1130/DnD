
cran_pkgs <- c("BiocManager", "ggplot2", "reshape2")
for (cran_pkg in cran_pkgs) {
  if (!require(cran_pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
  }
}

bioc_pkgs <- c("rtracklayer")
for (bioc_pkg in bioc_pkgs) {
  if (!require(bioc_pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(bioc_pkg)
  }
}

