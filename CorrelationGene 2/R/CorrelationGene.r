#' @title CorrelationGene
#'
#' @description Functions for CorrelationGene
#'
#' @param obj Seurat object
#' @param gene_list List of genes
#' @param genecolumn Column name of gene list
#' @param goi Gene of interest
#' @param celltype Celltype of interest
#'
#' @return Data frame with correlation values
#'
#' @examples
#' present_genes_RNA(obj, gene_list, genecolumn)
#' present_genes_SCT(obj, gene_list, genecolumn)
#' average_expression_SCT(obj, goi, gene_list)
#' average_expression_RNA(obj, goi, gene_list)
#' celltype_expression_SCT(obj, celltype, GOI, gene_list)
#' celltype_expression_RNA(obj, celltype, GOI, gene_list)
#' correlation(df, goi)
#' express_cell_RNA(obj, goi, gene_list)
#' express_cell_SCT(obj, goi, gene_list)
#' express_average_RNA(obj, goi, gene_list)
#' express_average_SCT(obj, goi, gene_list)
#' 
#' @export


library(Seurat)
library(glmGamPoi)
library(Matrix)
library(tictoc)
library(Hmisc)
library(moments)


SFARI_genes.names = load(file = "/Users/artemis/Sweden/CorrelationGene0.1.0/data/SFARI_genes.rda")


present_genes_RNA = function(obj, gene_list, genecolumn) {
  present_genes <- intersect(gene_list$genecolumn, rownames(obj@assays$RNA@data))
  return(present_genes)
}

present_genes_SCT = function(obj, gene_list, genecolumn) {
  present_genes <- intersect(gene_list$genecolumn, rownames(obj@assays$SCT@data))
  return(present_genes)
}

average_expression_SCT <- function(obj, goi, gene_list) {
    average_expression = AverageExpression(obj, features = c(gene_list,goi), layer = "data")$SCT
    average_expression_df <- as.data.frame(average_expression)
    average_expression_df = t(average_expression_df)
    View(average_expression_df)
    return(average_expression_df)
}

average_expression_RNA <- function(obj, goi, gene_list) {
    average_expression = AverageExpression(obj, features = c(gene_list,goi), layer = "data")$RNA
    average_expression_df <- as.data.frame(average_expression)
    average_expression_df = t(average_expression_df)
    View(average_expression_df)
    return(average_expression_df)
}


celltype_expression_SCT = function(obj, celltype, GOI, gene_list) {
    celltype_expression = subset(obj, idents = celltype)@assays$SCT@data
    celltype_expression = as.data.frame(celltype_expression)
    celltype_expression$gene = rownames(celltype_expression)
    celltype_expression_xbox <- celltype_expression[celltype_expression$gene %in% c(gene_list, GOI), ]
    celltype_expression_xbox = t(celltype_expression_xbox)
    celltype_expression_xbox = celltype_expression_xbox[-which(rownames(celltype_expression_xbox) == "gene"), ]
    return(celltype_expression_xbox)
}

celltype_expression_RNA = function(obj, celltype, GOI, gene_list) {
    celltype_expression = subset(obj, idents = celltype)@assays$RNA@data
    celltype_expression = as.data.frame(celltype_expression)
    celltype_expression$gene = rownames(celltype_expression)
    celltype_expression_xbox <- celltype_expression[celltype_expression$gene %in% c(gene_list, GOI), ]
    celltype_expression_xbox = t(celltype_expression_xbox)
    celltype_expression_xbox = celltype_expression_xbox[-which(rownames(celltype_expression_xbox) == "gene"), ]
    return(celltype_expression_xbox)
}

correlation = function(df, goi, gene_list) {
  cat("Computing correlation matrix...\n")

  df = rcorr(df, type = c("pearson"))
  df_r = as.data.frame(df$r)
  df_p = as.data.frame(df$P)
  cat("extracting data...\n")
  goi_only_p = as.data.frame(df_p[goi], row.names = rownames(df_p))
  goi_only_r = as.data.frame(df_r[goi], row.names = rownames(df_r))
  cat("adding gene names...\n")
  goi_only_p$gene = rownames(goi_only_p)
  goi_only_r$gene = rownames(goi_only_r)
  merged = merge(goi_only_r, goi_only_p, by = "gene")
  colnames(merged)[2:3] <- c("R","P")  
  merged = merged[order(merged$R, decreasing = TRUE), ]
  merged$SFARI.Gene = ""
  cat("Merging data frames...\n")
  common_genes = intersect(merged$gene, gene_list)
  
  is_sfari <- ifelse(SFARI_genes.names %in% common_genes, TRUE, FALSE)
  is_sfari <- is_sfari[1:nrow(merged)]

  cat("Adding gene annotations...\n")
  for (i in 1:nrow(merged)) {
    if (!is.na(is_sfari[i]) && is_sfari[i]) {
      merged[i, "SFARI.Gene"] <- "Y"
    }
  }
  cat("Done.\n")
  return(merged)
}

express_cell_RNA = function(obj, goi, gene_list) {
    merged = data.frame()
    tic()
    for (celltype in levels(obj)) {
        cat("Processing celltype: ", celltype, "into a dataframe", "\n")
        celltype_df_a = celltype_expression_RNA(obj, celltype, goi, gene_list)
        celltype_df_a = correlation(celltype_df_a, goi, gene_list)
        cat("Removing ", goi, " from celltype: ", celltype, "\n")
        celltype_a = celltype_df_a[!grepl(goi, celltype_df_a$gene),]
        cat("Adding celltype column...\n")
        celltype_a$celltype = celltype
        cat("Merging data frames...\n")
        merged = rbind(merged, celltype_a)
    }
    merged = merged[order(merged$R, decreasing = TRUE), ]
    cat("Done.\n")
    toc()

    
    return(merged)
}
express_cell_SCT = function(obj, goi, gene_list) {
    merged = data.frame()
    tic()
    for (celltype in levels(obj)) {
        cat("Processing celltype: ", celltype, " into a dataframe", "\n")
        celltype_df_a = celltype_expression_SCT(obj, celltype, goi, gene_list)
        celltype_df_a = correlation(celltype_df_a, goi, gene_list)
        cat("Removing ", goi, " from celltype: ", celltype, "\n")
        celltype_a = celltype_df_a[!grepl(goi, celltype_df_a$gene),]
        cat("Adding celltype column...\n")
        celltype_a$celltype = celltype
        cat("Merging data frames...\n")
        merged = rbind(merged, celltype_a)
    }
    merged = merged[order(merged$R, decreasing = TRUE), ]
    cat("Done.\n")
    toc()

    
    return(merged)
}
express_average_RNA = function(obj, goi, gene_list) {
    tic()
    cat("Processing gene: ", goi, "\n")
    average = average_expression_RNA(obj, goi, gene_list)
    cat("Computing correlation matrix...\n")
    average = correlation(average, goi)
    cat("Done.\n")
    toc()
    return(average)

}
express_average_SCT = function(obj, goi, gene_list) {
    tic()
    cat("Processing gene: ", goi, "\n")
    average = average_expression_SCT(obj, goi, gene_list)
    cat("Computing correlation matrix...\n")
    average = correlation(average, goi)
    cat("Done.\n")
    toc()
    return(average)

}





