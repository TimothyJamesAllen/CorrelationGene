#' Finds correlation between a gene of interest and all other genes in a Seurat object per level
#'
#' @param obj Seurat object. gene_list Gene list. goi Gene of interest. ensembl boolean
#' 
#' @return data frame with correlation values
#'
#' @examples
#' express_cell_RNA(obj, goi, gene_list)
#' 
#' @export
express_cell_RNA = function(obj, goi, gene_list, PorS, ensembl) {
    merged = data.frame()
  tic()
  for (celltype in levels(obj)) {
    cat("Processing celltype: ", celltype, " into a dataframe", "\n")
    celltype_df_a = celltype_expression_RNA(obj, celltype, goi, gene_list)
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


equation = function(a,x,R, padj) {
  b = a/1e+9
  c = a/2e+4
  y = a / (1+exp(-(b*x-c)))
  score = (y*R)/(log(padj+2))
  score = orderNorm(score)
  score = score$x.t
  return(score)
}
  cat = "Calculating scores...\n"

  a = max(merged$cell_number, na.rm = TRUE)
  print(a)

  merged$score = equation(a, merged$cell_number, merged$R, merged$padj)
  
  merged = merged[order(merged$score, decreasing = TRUE), ]



  toc()

  return(merged)
}

