#' Transfers labels from reference data to query data
#' 
#' Generates predictions for query data from reference data labels by transferring labels from reference to query data
#' Parameters:
#' @param query     Seurat object to query for predictions
#' @param ref       Seurat object to derive predictions from
#' @param refdata   Reference labels to use for prediction (cell type, cluster identity, other metadata, etc)
#' 
#' @import  Seurat
#' @return  Outputs query Seurat object with predictions added as metadata
#' @export

label_transfer <- function(query, ref, refdata){
    anchors = Seurat::FindTransferAnchors(reference = ref, query = query, dims = 1:30, project.query = T)
    predictions = Seurat::TransferData(anchorset = anchors, refdata = refdata, dims = 1:30)
    query = Seurat::AddMetaData(query, metadata = predictions)
    return(query)
}