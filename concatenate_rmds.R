# Script to concatenate all 13 phase Rmd files into one combined Rmd
# Run: Rscript concatenate_rmds.R

phases <- c(
  "phase1_scRNA_QC_preprocessing.Rmd",
  "phase2_normalization_dimreduction.Rmd",
  "phase3_clustering_annotation.Rmd",
  "phase4_monocyte_deepdive.Rmd",
  "phase5_replication_GSE150861.Rmd",
  "phase6_scRNAseq_visualization_summary.Rmd",
  "phase7_bulk_loading_EDA.Rmd",
  "phase8_DESeq2_differential_expression.Rmd",
  "phase9_adjust_confounders.Rmd",
  "phase10_bulk_visualization.Rmd",
  "phase11_PCA_clustering_bulk.Rmd",
  "phase12_GSEA_pathway_analysis.Rmd",
  "phase13_monocyte_signature_bulk.Rmd"
)

prefixes <- paste0("p", 1:13)

# Master YAML + setup
header <- '---
title: "COVID-19 scRNA-seq & Bulk RNA-seq Analysis - All Phases Combined"
author: ""
date: "`r Sys.Date()`"
output: html_document
---

# Combined Analysis Pipeline

This document contains **all 13 phases** of the COVID-19 PBMC analysis,
combined into a single executable notebook.

---

## Master Setup

```{r setup, message=FALSE, warning=FALSE}
library(Seurat)
library(Matrix)
library(tidyverse)
library(patchwork)
library(SingleR)
library(celldex)
library(dittoSeq)
library(pheatmap)
library(EnhancedVolcano)
library(ggrepel)
library(RColorBrewer)
library(DESeq2)
library(GEOquery)
library(org.Hs.eg.db)
library(apeglm)
library(factoextra)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(GSVA)

knitr::opts_chunk$set(
  echo    = TRUE,
  warning = FALSE,
  message = FALSE,
  fig.width  = 12,
  fig.height = 7
)
```

'

out_lines <- strsplit(header, "\n")[[1]]

for (i in seq_along(phases)) {
  fname <- phases[i]
  prefix <- prefixes[i]

  if (!file.exists(fname)) {
    cat("WARNING: File not found:", fname, "\n")
    next
  }

  lines <- readLines(fname, warn = FALSE)

  # Remove YAML front matter (everything between first --- and second ---)
  yaml_start <- which(lines == "---")
  if (length(yaml_start) >= 2) {
    yaml_end <- yaml_start[2]
    lines <- lines[(yaml_end + 1):length(lines)]
  }

  # Process setup chunk: remove library() calls and knitr::opts_chunk$set()
  # but keep other code (dir.create, variable definitions, etc.)
  in_setup <- FALSE
  setup_start <- NULL
  setup_end <- NULL
  for (j in seq_along(lines)) {
    if (!in_setup && grepl("^```\\{r\\s+setup", lines[j])) {
      in_setup <- TRUE
      setup_start <- j
    } else if (in_setup && grepl("^```$", lines[j])) {
      setup_end <- j
      break
    }
  }

  if (!is.null(setup_start) && !is.null(setup_end)) {
    # Extract the setup chunk body (between header and closing ```)
    setup_body <- lines[(setup_start + 1):(setup_end - 1)]
    
    # Keep lines that are NOT library(), knitr::opts_chunk$set, or their continuations
    keep_lines <- character(0)
    in_knitr_block <- FALSE
    for (sl in setup_body) {
      if (grepl("^\\s*knitr::opts_chunk", sl)) {
        in_knitr_block <- TRUE
        next
      }
      if (in_knitr_block) {
        if (grepl("^\\s*\\)", sl)) {
          in_knitr_block <- FALSE
          next
        }
        next
      }
      if (grepl("^\\s*library\\(", sl)) next
      if (grepl("^\\s*#.*", sl) && length(keep_lines) == 0) next  # skip leading comments
      if (grepl("^\\s*$", sl) && length(keep_lines) == 0) next    # skip leading blank lines
      keep_lines <- c(keep_lines, sl)
    }
    
    # Remove trailing blank lines
    while (length(keep_lines) > 0 && grepl("^\\s*$", keep_lines[length(keep_lines)])) {
      keep_lines <- keep_lines[-length(keep_lines)]
    }
    
    # Replace the setup chunk: if there's leftover code, make it a prefixed init chunk
    if (length(keep_lines) > 0) {
      init_chunk <- c(
        paste0("```{r ", prefix, "-init, message=FALSE, warning=FALSE}"),
        keep_lines,
        "```"
      )
      lines <- c(lines[1:(setup_start - 1)], init_chunk, lines[(setup_end + 1):length(lines)])
    } else {
      lines <- lines[-(setup_start:setup_end)]
    }
  }

  # Prefix all chunk labels to make them unique
  # Match ```{r chunk-name, ...} or ```{r chunk-name}
  lines <- gsub(
    "^```\\{r\\s+([^,}]+)",
    paste0("```{r ", prefix, "-\\1"),
    lines
  )

  # Add phase separator
  phase_title <- sub("_", " - ", sub("\\.Rmd$", "", fname))
  phase_title <- gsub("_", " ", phase_title)

  separator <- c(
    "",
    paste0("---"),
    "",
    paste0("# ===================================================================="),
    paste0("# ", toupper(phase_title)),
    paste0("# ===================================================================="),
    ""
  )

  out_lines <- c(out_lines, separator, lines)

  cat("Processed:", fname, "->", length(lines), "lines\n")
}

# Add final session info
out_lines <- c(out_lines, "", "---", "",
               "# Final Session Info", "",
               "```{r final-session-info}",
               "sessionInfo()",
               "```")

# Write output
writeLines(out_lines, "all_phases_combined.Rmd")
cat("\nDone! Written", length(out_lines), "lines to all_phases_combined.Rmd\n")
