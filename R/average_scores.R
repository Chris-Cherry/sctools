#' Averages values from a matrix based on selected metadata from Seurat object.
#' Matrix input should be output from organize_scores.
#'
#' @param scores    Matrix of scores organized by metadata (output from organize_scores)
#' @param meta      Selected metadata from Seurat object over which to calculate averages
#'
#' @return Matrix of average values 
#' @export

average_scores <- function(scores, meta){
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
  return(org_ave)
}