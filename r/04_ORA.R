rm(list = ls())
load(file = 'step1output.Rdata')

library(rtracklayer)
library(dplyr)
library(data.table)
library(fgsea)
library(ggplot2)

gtf_file <- "data/gencode.v46.chr_patch_hapl_scaff.annotation.gtf"
gtf_data <- import(gtf_file)

gtf_df <- as.data.frame(gtf_data)

gene_info <- gtf_df %>%
  filter(type == "gene") %>%
  dplyr::select(gene_id = gene_id, gene_name = gene_name)

gene_info$gene_id <- sub("\\..*", "", gene_info$gene_id)

important_features <- read.csv("data/feature_importance.csv")

important_features <- important_features %>%
  left_join(gene_info, by = c("feature" = "gene_id"))

read_gmt <- function(file) {
  gmt_lines <- readLines(file)
  gmt_list <- list()
  
  for (line in gmt_lines) {
    parts <- unlist(strsplit(line, "\t")) 
    gene_set_name <- parts[1] 
    genes <- parts[3:length(parts)] 
    
    genes <- genes[genes != ""]
    
    gmt_list[[gene_set_name]] <- genes 
  }
  
  return(gmt_list)
}

gmt_file <- "data/combined_ensembl.gmt"
gmt_data <- read_gmt(gmt_file)

res <- DEGs
ranks <- res$log2FoldChange
names(ranks) <- rownames(res)

ranks <- sort(ranks, decreasing = TRUE)
set.seed(17)
fgsea_res <- fgsea(pathways = gmt_data, 
                   stats    = ranks,
                   eps      = 0.0,
                   minSize  = 5,
                   maxSize  = 500)

fgsea_res <- fgsea_res[fgsea_res$padj < 0.01, ]


top_genes <- important_features$feature[1:20]

ora_test <- function(gene_set, genes_of_interest, background_genes) {
  
  overlap <- length(intersect(gene_set, genes_of_interest))
  set_size <- length(gene_set)
  total_genes <- length(background_genes)
  target_size <- length(genes_of_interest)
  
  p_value <- phyper(overlap - 1, set_size, total_genes - set_size, target_size, lower.tail = FALSE)
  return(p_value)
}

background_genes <- rownames(DEGs)

ora_results <- data.frame(GeneSet = character(), PValue = numeric(), stringsAsFactors = FALSE)

for (gene_set_name in names(gmt_data)) {
  gene_set <- gmt_data[[gene_set_name]]
  p_value <- ora_test(gene_set, top_genes, background_genes)
  
  ora_results <- rbind(ora_results, data.frame(GeneSet = gene_set_name, PValue = p_value))
}

ora_results$AdjustedPValue <- p.adjust(ora_results$PValue, method = "BH")

ora_results$GeneSetSize <- sapply(gmt_data[ora_results$GeneSet], length)
ora_results$Overlap <- sapply(ora_results$GeneSet, function(gs) length(intersect(gmt_data[[gs]], top_genes)))
ora_results$EnrichmentFactor <- ora_results$Overlap / (ora_results$GeneSetSize * length(top_genes) / length(background_genes))

ora_results$OverlapGenes <- sapply(ora_results$GeneSet, function(gs) {
  overlap_ids <- intersect(gmt_data[[gs]], top_genes)
  
  overlap_names <- gene_info %>%
    filter(gene_id %in% overlap_ids) %>%
    pull(gene_name)
  
  paste(overlap_names, collapse = ", ")
})

ora_results <- ora_results %>% arrange(AdjustedPValue)

significant_results <- ora_results %>% filter(AdjustedPValue < 0.01)
print(significant_results)

significant_results$source <- ifelse(grepl("R-HSA-|Homo sapiens", significant_results$GeneSet), "Reactome",
                           ifelse(grepl("WP", significant_results$GeneSet), "WikiPathways",
                                  ifelse(grepl("_|outerRadialGlia", significant_results$GeneSet), "Literature",
                                         ifelse(grepl("UP|DN|\\.V[1-9]|CAHOY|PID PLK1 PATHWAY", significant_results$GeneSet) & !grepl("DNA", significant_results$GeneSet), "MsigDB",
                                                "KEGG"))))

significant_results <- significant_results %>%
  mutate(GeneSet = case_when(
    GeneSet == "Glial Cell Differentiation WP2276" ~ "Glial Cell Differentiation (WP2276)",
    TRUE ~ GeneSet
  ))
