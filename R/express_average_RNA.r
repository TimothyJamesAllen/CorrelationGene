#' Finds correlation between a gene of interest and all other genes across every level in a Seurat Obj
#'
#' 
#'
#' @param obj Seurat object. gene_list Gene list. goi Gene of interest.
#' 
#' @return data frame with correlation values
#'
#' @examples
#' express_average_RNA(obj, goi, gene_list)
#' 
#'
#' @export
express_average_RNA = function(obj, goi, gene_list) {
    tic()
    cat("Processing gene: ", goi, "\n")
    average = average_expression_RNA(obj, goi, gene_list)
    cat("Computing correlation matrix...\n")
    average = correlation(average, goi)
    cat("Done.\n")
    toc()
    return(average)

}