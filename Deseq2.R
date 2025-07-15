library(openxlsx)
library(DESeq2)

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
  res <- res[!is.na(res$padj) & res$padj < 0.05, ]
  
  fc_1.5 <- data.frame()
  fc_2   <- data.frame()
  fc_4   <- data.frame()
  
  for (gene in rownames(res)) {
    untreated_vals <- norm_counts[gene, replicates[[untreated]]]
    treated_vals   <- norm_counts[gene, replicates[[treated]]]
    
    if (max(untreated_vals) < 5) next
    
    fold_changes <- treated_vals / untreated_vals
    
    gene_info <- data.frame(
      Gene = gene,
      padj = res[gene, "padj"],
      FoldChanges = paste0(round(fold_changes, 2), collapse = ", ")
    )
    
    if (all(fold_changes > 1.5)) fc_1.5 <- rbind(fc_1.5, gene_info)
    if (all(fold_changes > 2))   fc_2   <- rbind(fc_2, gene_info)
    if (all(fold_changes > 4))   fc_4   <- rbind(fc_4, gene_info)
  }
  
  return(list(fc_1.5 = fc_1.5, fc_2 = fc_2, fc_4 = fc_4))
}

treated_m1    <- grep("^treated.*_m1$", names(replicate_groups), value = TRUE)
untreated_m1  <- grep("^untreated.*_m1$", names(replicate_groups), value = TRUE)
treated_m10   <- grep("^treated.*_m10$", names(replicate_groups), value = TRUE)
untreated_m10 <- grep("^untreated.*_m10$", names(replicate_groups), value = TRUE)

fc1.5_list <- list()
fc2_list   <- list()
fc4_list   <- list()

for (treated in c(treated_m1, treated_m10)) {
  untreated_group <- if (grepl("_m1$", treated)) untreated_m1 else untreated_m10
  for (untreated in untreated_group) {
    treated_name <- gsub("^treated_", "", treated)
    untreated_name <- gsub("^untreated_", "", untreated)
    name <- paste(treated_name, "vs", untreated_name, sep = "_")
    cat("Running:", name, "\n")
    
    fc_tables <- get_fc_tables(treated, untreated, norm_counts, dds, replicate_groups)
    
    if (nrow(fc_tables$fc_1.5) > 0) fc1.5_list[[name]] <- fc_tables$fc_1.5
    if (nrow(fc_tables$fc_2) > 0)   fc2_list[[name]]   <- fc_tables$fc_2
    if (nrow(fc_tables$fc_4) > 0)   fc4_list[[name]]   <- fc_tables$fc_4
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

write_fc_workbook(fc1.5_list, "foldchange_1.5.xlsx")
write_fc_workbook(fc2_list,   "foldchange_2.xlsx")
write_fc_workbook(fc4_list,   "foldchange_4.xlsx")


