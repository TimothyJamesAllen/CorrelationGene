#' Finds correlation between a gene of interest and all other genes in a Seurat object
#'
#' 
#'
#' 
#' @param df Matrix of gene expression values
#' @param PorS Type of correlation to compute (P for Pearson, S for Spearman)
#' @param goi Gene of interest.
#' @param gene_list Gene list.
#' 
#' 
#' @return A matrix of correlation values
#'
#' @examples
#' correlation(df, goi, gene_list, PorS)
#' correlation(matrix,"RFX3", gene_list, "P")
#' 
#'
#' @export
correlation = function(df, goi, gene_list, PorS) {

  cat("Computing correlation matrix...\n")
  if (PorS == "P") {
    df <- rcorr(as.matrix(df), type = "pearson")
  } else if (PorS == "S") {
    df <- rcorr(as.matrix(df), type = "spearman")
  } else {
    cat("Error: PorS must be either P or S")
  }
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
  
  is_sfari <- ifelse(SFARI_genes$gene.symbol %in% common_genes, TRUE, FALSE)
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

