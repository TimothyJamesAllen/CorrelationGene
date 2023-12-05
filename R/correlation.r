#' Finds correlation between a gene of interest and all other genes in a Seurat object
#'
#' 
#'
#' @param obj Seurat object. gene_list Gene list. goi Gene of interest.
#' 
#' @return A vector of present genes
#'
#' @examples
#' correlation(df, goi, gene_list)
#' 
#'
#' @export
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