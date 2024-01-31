#' Creates expression matrix for a gene of interest and a list of genes and a celltype
#'
#' 
#'
#' @param obj Seurat object. goi Gene of interest. gene_list Gene list. celltype Celltype of interest
#'
#' @return A vector of present genes
#'
#' @examples
#' celltype_expression_RNA(obj, celltype, GOI, gene_list)
#' 
#'
#' @export
celltype_expression_RNA = function(obj, celltype, GOI, gene_list) {
    celltype_expression = subset(obj, idents = celltype)@assays$RNA@data
    celltype_expression = as.data.frame(celltype_expression)
    celltype_expression$gene = rownames(celltype_expression)
    celltype_expression_xbox <- celltype_expression[celltype_expression$gene %in% c(gene_list, GOI), ]
    celltype_expression_xbox = t(celltype_expression_xbox)
    celltype_expression_xbox = celltype_expression_xbox[-which(rownames(celltype_expression_xbox) == "gene"), ]
    return(celltype_expression_xbox)
}

