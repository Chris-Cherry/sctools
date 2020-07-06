#' Sort matrix by Seurat metadata
#' 
#' Organizes matrix of values associated with cells in a Seurat object 
#' by specified metadata
#'
#' @param scores    The matrix of values to be organized
#' @param ser       The Seurat object containing the cells of interest
#' @param meta      The metadata used to organize the values
#'
#' @return Organized matrix of scores
#' @export

organize_scores <- function(scores, ser, meta){
  org = meta
  names(org) = colnames(ser@assays$RNA@counts)
  org = sort(org)
  if (is.null(nrow(scores))){
    org_score = scores[match(names(org), names(scores))]
  } else {
    org_score = scores[,match(names(org), colnames(scores))]
  }
  return(org_score)
}
