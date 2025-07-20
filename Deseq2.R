# ORFeome – Analyzing ORFeome screening data
# Copyright (C) 2025 Anushka Shome and Franck Dumetz
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see https://www.gnu.org/licenses/.

# This script identifies differentially expressed genes from ORFeome screening data and exports fold-change results to Excel.

library(openxlsx)
library(DESeq2)

args <- commandArgs(trailingOnly = TRUE)
fc_input <- as.integer(args[1])

counts <- read.csv("counts.csv", header = TRUE, row.names = 1)
counts <- as.matrix(counts)

coldata <- read.csv("treatments.csv", row.names = 1, stringsAsFactors = FALSE)
coldata$Condition <- factor(coldata$Condition)

stopifnot(all(rownames(coldata) == colnames(counts)))

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~Condition)
dds <- DESeq(dds)
norm_counts <- counts(dds, normalized = TRUE)

replicate_groups <- split(colnames(norm_counts), coldata$Condition)

get_fc_tables <- function(treated, untreated, norm_counts, dds, replicates) {
  contrast <- c("Condition", treated, untreated)
  res <- results(dds, contrast = contrast)
  
  fc <- data.frame()
  
  res <- res[!is.na(res$padj) & res$padj < 0.05, ]
  
  for (gene in rownames(res)) {
    untreated_vals <- norm_counts[gene, replicates[[untreated]]]
    treated_vals   <- norm_counts[gene, replicates[[treated]]]
    
    fold_changes <- treated_vals / untreated_vals
    
    gene_info <- data.frame(
      Gene = gene,
      padj = res[gene, "padj"],
      FoldChanges = paste0(round(fold_changes, 2), collapse = ", ")
    )
    
    if (max(untreated_vals) < 5) next
    
    if (all(fold_changes > fc_input)) fc <- rbind(fc, gene_info)
  }
  
  return(list(fc=fc))
}

no_counts <- function(untreated, replicates){
  zero <- data.frame()
  untreated <- unlist(replicates[untreated])
  for (gene in rownames(norm_counts)) {
    untreated_vals <- norm_counts[gene, untreated]
    if (max(untreated_vals) < 5) zero <- rbind(zero, gene)
  }
  return(list(zero))
}

treated_m1    <- grep("^treated.*_m1$", names(replicate_groups), value = TRUE)
untreated_m1  <- grep("^untreated.*_m1$", names(replicate_groups), value = TRUE)
treated_m10   <- grep("^treated.*_m10$", names(replicate_groups), value = TRUE)
untreated_m10 <- grep("^untreated.*_m10$", names(replicate_groups), value = TRUE)
untreated_both <- grep("^untreated", names(replicate_groups), value = TRUE)

zero <- list()
fc_list <- list()

zero <- no_counts(untreated_both, replicate_groups)

for (treated in c(treated_m1, treated_m10)) {
  untreated_group <- if (grepl("_m1$", treated)) untreated_m1 else untreated_m10
  for (untreated in untreated_group) {
    treated_name <- gsub("^treated_", "", treated)
    untreated_name <- gsub("^untreated_", "", untreated)
    name <- paste(treated_name, "vs", untreated_name, sep = "_")
    cat("Running:", name, "\n")
    
    fc_tables <- get_fc_tables(treated, untreated, norm_counts, dds, replicate_groups)
    
    fc_list[[name]] <- fc_tables$fc
  }
}

write_fc_workbook <- function(fc_list, filename) {
  wb <- createWorkbook()
  
  if (length(fc_list) > 0) {
    common_genes <- Reduce(intersect, lapply(fc_list, function(x) x$Gene))
    
    if (length(common_genes) > 0) {
      common_df <- data.frame(Gene = common_genes)
      addWorksheet(wb, "CommonHits")
      writeData(wb, "CommonHits", common_df)
    }
    
    for (name in names(fc_list)) {
      addWorksheet(wb, name)
      writeData(wb, name, fc_list[[name]])
    }
  }
  
  saveWorkbook(wb, filename, overwrite = TRUE)
}

write_fc_workbook(fc_list, paste0("foldchange_", fc_input, ".xlsx"))

zero_wb <- createWorkbook()

addWorksheet(zero_wb, "no-counts")

writeData(zero_wb, "no-counts", zero)

saveWorkbook(zero_wb, "no-counts.xlsx", overwrite = TRUE)

zero_wb <- createWorkbook()
addWorksheet(zero_wb, "no-counts")
writeData(zero_wb, "no-counts", zero)
saveWorkbook(zero_wb, "no-counts.xlsx", overwrite = TRUE)
