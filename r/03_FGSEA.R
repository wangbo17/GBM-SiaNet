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

names(gmt_data)[names(gmt_data) == "PLK1 signaling events Homo sapiens e5e87977-6194-11e5-8ac5-06603eb7f303"] <- "PID PLK1 PATHWAY"

res <- DEGs
ranks <- res$log2FoldChange
names(ranks) <- rownames(res)

ranks <- sort(ranks, decreasing = TRUE)
set.seed(17)
fgsea_res <- fgsea(pathways = gmt_data, 
                  stats    = ranks,
                  eps      = 0.0,
                  minSize  = 15,
                  maxSize  = 500)

fgsea_res$source <- ifelse(grepl("R-HSA-|Homo sapiens", fgsea_res$pathway), "Reactome",
                              ifelse(grepl("WP", fgsea_res$pathway), "WikiPathways",
                                     ifelse(grepl("_|outerRadialGlia", fgsea_res$pathway), "Literature",
                                            ifelse(grepl("UP|DN|\\.V[1-9]|CAHOY|PID PLK1 PATHWAY", fgsea_res$pathway) & !grepl("DNA", fgsea_res$pathway), "MsigDB",
                                                   "KEGG"))))
fgsea_res <- fgsea_res[fgsea_res$padj < 0.01, ]

top20_features <- important_features %>%
  dplyr::slice(1:20) %>%
  dplyr::select(feature)

fgsea_res$top20_genes_in_pathway <- NA
fgsea_res$num_top20_genes_in_pathway <- 0

for (i in 1:nrow(fgsea_res)) {
  pathway_genes <- gmt_data[[fgsea_res$pathway[i]]]
  
  matching_genes <- top20_features$feature[top20_features$feature %in% pathway_genes]
  
  fgsea_res$top20_genes_in_pathway[i] <- paste(matching_genes, collapse = ", ")
  fgsea_res$num_top20_genes_in_pathway[i] <- length(matching_genes)
}

fgsea_res$top20_genes_in_pathway <- sapply(fgsea_res$top20_genes_in_pathway, function(genes) {
  if (genes != "") {
    gene_ids <- unlist(strsplit(genes, ", "))
    gene_names <- gene_info$gene_name[match(gene_ids, gene_info$gene_id)]
    return(paste(gene_names, collapse = ", "))
  } else {
    return(NA)
  }
})

split_data <- split(fgsea_res, fgsea_res$source)

kegg_df <- split_data[["KEGG"]]
msigdb_df <- split_data[["MsigDB"]]
reactome_df <- split_data[["Reactome"]]
wikipathways_df <- split_data[["WikiPathways"]]
literature_df <- split_data[["Literature"]]
