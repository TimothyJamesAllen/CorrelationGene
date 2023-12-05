#' Creates average expression matrix for a gene of interest and a list of genes 
#'
#' 
#'
#' @param obj Seurat object. gene_list Gene list. genecolumn Gene column name
#'
#' @return A vector of present genes
#'
#' @examples
#' average_expression_SCT(obj, goi, gene_list)
#' 
#'
#' @export
average_expression_SCT <- function(obj, goi, gene_list) {
    average_expression = AverageExpression(obj, features = c(gene_list,goi), layer = "data")$SCT
    average_expression_df <- as.data.frame(average_expression)
    average_expression_df = t(average_expression_df)
    View(average_expression_df)
    return(average_expression_df)
}