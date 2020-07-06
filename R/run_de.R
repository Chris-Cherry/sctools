#' Runs differential expression for a ser object
#'
#' Creates object used in make_de_plots function for the creation of various plots.
#' Can be used on whole data set or individual clusters.
#' @param ser           The Seurat object to use
#' @param feats         A subset of featurs to run DE (I.E. only surface markers)
#'                      If NULL then will run on all genes
#' @param out_dir       Directory to write plots and save markers.
#' @param meta          (Optional) Pass in cluster name that user wants to run DE
#' @param origin        (Optional) Pass in what parameters in the each clusters that user wants to run DE
#' @param verbose          Boolean flag for verbose output
#' @import grDevices
#' @import Seurat
#' @importFrom utils write.table
#' @importFrom dplyr %>% 
#' @return              Output a processed differential expression for plotting
#' @export

run_de <- function(ser, feats = NULL, out_dir = '2_de/', meta = NULL, origin = 'Condition', verbose = F){
    dir.create(out_dir)
    out_tmp = out_dir

    # Run DE on whole dataset based on the meta data
    if(is.null(feats)){
        DefaultAssay(ser) = "RNA"
        feats = rownames(ser)}
    ser_list = list()
    j = 1
    # Run DE on each cluster
    if (!is.null(meta)){
        for(clust in meta){
            out_dir = out_tmp
            subser = subset(ser, idents = clust)
            Idents(subser) = subser@meta.data[origin]
            out_dir = paste0(out_dir, '/cluster_', clust, '/')
            dir.create(out_dir)
            sub_marks = FindAllMarkers(subser, features = feats, logfc.thresh = 0,
                return.thresh = Inf, assay = 'RNA', verbose = verbose)
            saveRDS(sub_marks, paste0(out_dir, 'submarks_', clust,'.RDS'))
            ser_list[[clust]] <- subser
            j = j + 1

            for(clust in levels(Idents(subser))){
                cl_de = sub_marks[which(sub_marks$cluster == clust),]
                write.table(cl_de, paste0(out_dir, '/', clust, '_de.csv'), sep = ',', 
                    col.names = NA, quote = FALSE)
            }
        }
    return(ser_list)
    }

    else{
        marks = FindAllMarkers(ser, features = feats, logfc.thresh = 0, 
            return.thresh = Inf, assay = 'RNA', verbose = verbose)
        saveRDS(marks, paste0(out_dir, '/marks.RDS'))

        for(clust in levels(Idents(ser))){
        cl_de = marks[which(marks$cluster == clust),]
        write.table(cl_de, paste0(out_dir, '/', clust, '_de.csv'), sep = ',', 
            col.names = NA, quote = FALSE)
        }
    return(ser)
    }

    dev.off()
}
