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
  
  is_sfari <- rep(FALSE, nrow(merged))
  is_sfari[merged$gene %in% SFARI_genes$gene.symbol] <- TRUE

cat("Adding gene annotations...\n")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene info for all genes at once
gene_info <- getBM(attributes = c("external_gene_name", "description"), filters = "external_gene_name", values = merged$gene, mart = mart)

# Create a named vector of descriptions with gene names as names
descriptions <- setNames(gene_info$description, gene_info$external_gene_name)

# Add descriptions to merged data frame
merged$description <- ifelse(merged$gene %in% names(descriptions), descriptions[merged$gene], "N/A")

# Add "Y" to "SFARI.Gene" column for rows where is_sfari is TRUE
merged[is_sfari, "SFARI.Gene"] <- "Y"
  
cat("Removing GOI...\n")

merged <- merged[-which(merged$gene == goi), ]


cat("Done.\n")
  return(merged)
}

