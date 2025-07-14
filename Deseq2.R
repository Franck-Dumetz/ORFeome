#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")

library(DESeq2)

counts <- read.csv("counts.csv", header=TRUE, row.names=1)
counts <- as.matrix(counts)
#condition <- factor(c('NewMP_m1' , 'NewMP_m1' , 'NewMP_m1' , 'OldMP_m1','OldMP_m1' , 'OldMP_m1' ,'OneTX_m1' , 'OneTX_m1' , 'OneTX_m1' , 'TwoTX_m1' ,'TwoTX_m1', 'TwoTX_m1'))
#condition = factor(c("TwoTX_m1", "OneTX_m1", "OneTX_m1", "OneTX_m1", "NewMP_m1", "NewMP_m1", "OldMP_m1", "NewMP_m1", "TwoTX_m1", "TwoTX_m1", "OldMP_m1", "OldMP_m1"))
#coldata <- data.frame(condition=condition)

coldata <- read.csv("treatments.csv", row.names=1, stringsAsFactors = FALSE)
head(coldata)
coldata$Condition <- factor(coldata$Condition)

all(rownames(coldata) == colnames(counts))
diff(rownames(coldata) == colnames(counts))


#rownames(coldata) <- colnames(counts)
dim(counts)
dim(coldata)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~Condition)
dds <- DESeq(dds)

norm_counts <- counts(dds, normalized=TRUE)

replicate_groups <- split(colnames(norm_counts), coldata$Condition)

get_hits <- function(treated, untreated, norm_counts, dds, replicates) {
  contrast <- c("Condition", treated, untreated)
  res <- results(dds, contrast=contrast)
  sig_genes <- rownames(res)[!is.na(res$padj) & res$padj < 0.05]
  
  filtered_genes <- c()
  
  for (gene in sig_genes) {
    untreated_vals <- norm_counts[gene, replicates[[untreated]]]
    treated_vals <- norm_counts[gene, replicates[[treated]]]
    
    if (max(untreated_vals) < 5) next
    
    fold_changes <- treated_vals / untreated_vals
    if (all(fold_changes > 4)) {
      filtered_genes <- c(filtered_genes, gene)
    }
  }
  
  return(filtered_genes)
}

treated_m1 <- grep("^treated.*_m1$", names(replicate_groups), value = TRUE)
untreated_m1 <- grep("^untreated.*_m1$", names(replicate_groups), value = TRUE)
treated_m10 <- grep("^treated.*_m10$", names(replicate_groups), value = TRUE)
untreated_m10 <- grep("^untreated.*_m10$", names(replicate_groups), value = TRUE)

all_hits <- list()

for (treated in treated_m1) {
  for (untreated in untreated_m1) {
    hits <- get_hits(treated, untreated, norm_counts, dds, replicate_groups)
    name <- paste(treated, "vs", untreated, sep = "_")
    all_hits[[name]] <- hits
  }
}

for (treated in treated_m10) {
  for (untreated in untreated_m10) {
    hits <- get_hits(treated, untreated, norm_counts, dds, replicate_groups)
    name <- paste(treated, "vs", untreated, sep = "_")
    all_hits[[name]] <- hits
  }
}

str(all_hits)
final_hits <- Reduce(intersect, all_hits)
final_hits <- data.frame(Genes = final_hits)
write.csv(final_hits, file = "overrep_genes.csv", row.names = FALSE, quote = FALSE)
