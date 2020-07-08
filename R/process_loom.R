#' Reads in loom files for rnavel_plot.R
#' 
#' The function reads in loom files from align_rnavel.sh and formats the cell
#' names to match the naming convention of process_ser.R. The input files must
#' have the 'sample.loom' where sample is the correct sample name to append to
#' the cell barcode. Note that cells in the Seurat object must be named in the format
#' barcode_sample.
#' 
#' @param file Path to the output file of align_rnavel.sh. The name must match sample names in the provided metadata.csv file IF hashing was not performed. If hashing was performed then sample names will be determined from the hash meta data csv provided.
#' @param ser Seurat object containing target cells. If hashing is TRUE then it is recommended that ONLY the hashed cells are included in the given Seurat object. This is to prevent barcode overlap from different sequencing runs when searching for barcodes without the sample appended.
#' @param hashed Boolean indicating whether hashing was used to determine cell sample. If it was, then the hashing metadata file must be provided for the correct cell barcode format.
#' @return Outputs a list of loom files from RNA velocity prepared for rnavel_plot.R or concatenation if processing multiple samples.
#' @export
#'

process_loom <- function(file, ser, hashed = FALSE){
    loom = velocyto.R::read.loom.matrices(file)
    bcs = sapply(colnames(loom$spliced), function(x){
            x = strsplit(x, ':', fixed = TRUE)[[1]][2]
            return(substr(x, 1, nchar(x)-1))  
    })
    file_name = sapply(file, function(x){
            x = strsplit(x, '/', fixed = TRUE)
            x = x[[1]]
            x = x[length(x)]})
    file_name = sapply(file_name, function(x){
            x = strsplit(x, '.', fixed = TRUE)[[1]][1]})
    bcs = paste0(bcs,"_",file_name)
    names(bcs) = c()

    # Fixing the unconventional column name
    tar_bcs = colnames(ser@assays$RNA@counts)
    names(tar_bcs) = tar_bcs

    if(hashed){
        tar_bcs = sapply(tar_bcs, function(x){
            strsplit(x, '-', fixed = TRUE)[[1]][2]
        })
    }

    overlap = intersect(bcs, tar_bcs)
    loom = lapply(loom, function(x){
        colnames(x) = bcs
        x = x[,overlap]
        colnames(x) = names(tar_bcs)[match(overlap, tar_bcs)]
        return(x)
    })
    return(loom)
}
