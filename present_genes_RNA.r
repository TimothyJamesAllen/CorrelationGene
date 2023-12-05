#' Find Present genes in the seurat object and the gene list 
#'
#' 
#'
#' @param obj Seurat object. gene_list Gene list. genecolumn Gene column name
#'
#' @return A vector of present genes
#'
#' @examples
#' present_genes_RNA(obj, gene_list, genecolumn)
*' 
#'
#' @export
present_genes_RNA = function(obj, gene_list, genecolumn) {
  present_genes <- intersect(gene_list$genecolumn, rownames(obj@assays$RNA@data))
  return(present_genes)
}

