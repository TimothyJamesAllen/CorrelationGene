#' Finds correlation between a gene of interest and all other genes in a Seurat object per level
#'
#' 
#'
#' @param obj Seurat object. gene_list Gene list. goi Gene of interest.
#' 
#' @return data frame with correlation values
#'
#' @examples
#' express_cell_RNA(obj, goi, gene_list)
#' 
#'
#' @export
express_cell_RNA = function(obj, goi, gene_list, PorS) {
    merged = data.frame()
    tic()
    for (celltype in levels(obj)) {
        cat("Processing celltype: ", celltype, "into a dataframe", "\n")
        celltype_df_a = celltype_expression_RNA(obj, celltype, goi, gene_list)
        if (PorS == "P") {
            celltype_df_a <- correlation(celltype_df_a, goi, gene_list, "P")
        } else if (PorS == "S") {
            celltype_df_a <- correlation(celltype_df_a, goi, gene_list, "S")
        } else {
            cat("Error: PorS must be either P or S")
        }
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