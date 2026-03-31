# Ensure user library exists
user_lib <- Sys.getenv("R_LIBS_USER")
dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(user_lib)

# Fix download method for Windows
options(download.file.method = "wininet")

# CRAN packages
install.packages(c(
  "Seurat",
  "tidyverse",
  "patchwork",
  "Matrix",
  "BiocManager",
  "pheatmap",
  "ggrepel",
  "factoextra",
  "msigdbr",
  "devtools"
), repos = "https://cran.rstudio.com/", lib = user_lib)

# Bioconductor packages
BiocManager::install(c(
  "SingleR",
  "celldex",
  "SingleCellExperiment",
  "DESeq2",
  "apeglm",
  "org.Hs.eg.db",
  "clusterProfiler",
  "enrichplot",
  "GSVA",
  "dittoSeq",
  "EnhancedVolcano",
  "glmGamPoi"
), ask = FALSE, update = FALSE, lib = user_lib)
