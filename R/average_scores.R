#'Average matrix by metadata
#' 
#'Averages values from a matrix based on selected metadata from Seurat object.
#'
#' @param scores    Matrix of scores organized by metadata (output from organize_scores)
#' @param ser       Seurat object containing metadata and cell names
#' @param meta      Selected metadata from Seurat object over which to calculate averages
#'
#' @return Matrix of average values 
#' @export

average_scores <- function(scores, ser, meta){
  names(meta) = colnames(ser@assays$RNA@counts)
  if (is.null(nrow(scores))){
    org_ave = rep(NA, length(levels(meta)))
    names(org_ave) = levels(meta)
    for (i in 1:length(levels(meta))){
      org_ave[i] = mean(scores[match(names(meta[which(meta == levels(meta)[i])]), 
          names(scores))])
    }
  } else {
    org_ave = matrix(rep(NA, 
      length(rownames(scores))*length(levels(meta))), 
      nrow = length(levels(meta)))
    colnames(org_ave) = rownames(scores)
    rownames(org_ave) = levels(meta)
    for (i in 1:length(levels(meta))){
      for (j in 1:length(rownames(scores))){
        org_ave[i,j] = mean(scores[j, 
          match(names(meta[which(meta == levels(meta)[i])]), 
          colnames(scores))])
      }
    }
  }
  return(org_ave)
}
