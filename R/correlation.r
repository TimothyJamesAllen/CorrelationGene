#' Finds correlation between a gene of interest and all other genes in a Seurat object
#'
#' 
#'
#' 
#' @param df Matrix of gene expression values
#' @param PorS Type of correlation to compute (P for Pearson, S for Spearman)
#' @param goi Gene of interest.
#' @param gene_list Gene list.
#' @param ensembl Boolean, whether to convert Ensembl IDs to gene symbols 
#' 
#' 
#' @return A matrix of correlation values
#'
#' @examples
#' correlation(df, goi, gene_list, PorS, ensembl)
#' correlation(matrix,"RFX3", gene_list, "P", ensembl = "FALSE")
#' 
#'
#' @export
correlation = function(df, goi, gene_list, PorS, ensembl) {

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
  cat("Merging data frames...\n")
  common_genes = intersect(merged$gene, gene_list)
  merged = merged[merged$gene %in% common_genes, ]



  if (ensembl == TRUE) {
    cat("converting ensemble IDs to gene symbols...\n")
    symbols <- AnnotationDbi::select(org.Hs.eg.db, keys = merged$gene,
                                 columns = "SYMBOL",
                                 keytype = "ENSEMBL")

  }

  merged <- merge(merged, symbols, by.x = "gene", by.y = "ENSEMBL")

# Remove the original Ensembl gene ID column
  merged <- merged[, !(names(merged) %in% "gene")]

# Rename the new gene symbol column
  colnames(merged)[colnames(merged) == "SYMBOL"] <- "gene_symbol"


  cat("Done.\n")


# Loop through each gene ID in the merged data frame

for (i in 1:nrow(merged)) {
  n <- merged$external_gene_name[i]
  if (n %in% SFARI_genes$gene.symbol) {
    merged$SFARI.Gene[i] <- "TRUE"
  } else {
    merged$SFARI.Gene[i] <- "FALSE"
  }
}

cat("Removing", goi, "...\n")

merged <- merged[-which(merged$gene == goi), ]

cat("Done.\n")

return(merged)
}

