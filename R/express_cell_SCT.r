#' Finds correlation between a gene of interest and all other genes in a Seurat object per level
#'
#'
#'
#' @param obj Seurat object. gene_list Gene list. goi Gene of interest. ensembl boolean
#'
#' @return matrix
#'
#' @examples
#' express_cell_SCT(obj, goi, gene_list, PorS, ensembl)
#'
#'
#' @export
express_cell_SCT = function(obj, goi, gene_list, PorS, ensembl) {
    merged = data.frame()
  tic()
  for (celltype in levels(obj)) {
    cat("Processing celltype: ", celltype, " into a dataframe", "\n")
    celltype_df_a = celltype_expression_SCT(obj, celltype, goi, gene_list)
    print(dim(celltype_df_a))
    cell_number = dim(celltype_df_a)[1]
    if (is.null(cell_number) || cell_number == "") {
      cat("Skipping celltype: ", celltype, " because it has no cells", "\n")
      next
    }

  

    if (PorS == "P") {
      celltype_df_a <- correlation(celltype_df_a, goi, gene_list, "P", ensembl)
    } else if (PorS == "S") {
      celltype_df_a <- correlation(celltype_df_a, goi, gene_list, "S", ensembl)
    } else {
      cat("Error: PorS must be either P or S")
    }
    
    cat("Adding celltype column...\n")
    celltype_df_a$celltype = celltype
    celltype_df_a$cell_number = cell_number
    cat("Merging data frames...\n")
    merged = rbind(merged, celltype_df_a)
  }
  
  merged = merged[order(merged$R, decreasing = TRUE), ]
  cat("Calculated padj values...\n")
  merged$padj = p.adjust(merged$P, method = "BH")
  cat("Done.\n")


equation = function(x,R, padj) {
  a = 30000
  b = 6.9
  c = 0.000675
  y = a / (1+exp(-(b*x-c)))
  score = (y*R)/(log(padj+2))
  return(score)
}
  cat = "Calculating scores...\n"

  merged$score = equation(a, merged$cell_number, merged$R, merged$padj)
  
  merged = merged[order(merged$score, decreasing = TRUE), ]



  toc()

  return(merged)
}
