#' Runs differential expression and creates relevant plots for a ser object
#' Parameters:
#' @param ser           The Seurat object to use
#' @param feats         A subset of featurs to run DE (I.E. only surface markers)
#'                      If NULL then will run on all genes
#' @param out_dir       Directory to write plots and save markers.
#' @param meta          (Optional) Pass in meta data that user wants to subset
#' @param res           (Optional) Resolution for clustering in each subset
#' @import grDevices
#' @import Seurat
#' @import dplyr
#' @importFrom dplyr %>% 
#' @return              Output a processed differenttial expression for plotting
#' @export

run_de <- function(ser, feats = NULL, out_dir = '2_de/', meta = NULL, res = .8){
    dir.create(out_dir)

    # Run DE on whole dataset based on the meta data

    if(is.null(feats)){feats = rownames(ser)}

    ser_list = list()

    # Run DE on each cluster

    if (!is.null(meta)){
        for(clust in meta){
            subser = subset(ser, idents = clust)
            subser = FindNeighbors(subser, reduction = "mnn", verbose = FALSE)
            subser = FindClusters(subser, resolution = res, verbose = FALSE)
            out_dir = paste0('2_de/cluster_', clust, '/')
            dir.create(out_dir)
            sub_marks = FindAllMarkers(subser, features = feats, logfc.thresh = 0,
                return.thresh = Inf, assay = 'RNA', verbose = FALSE)
                saveRDS(sub_marks, paste0(out_dir, 'submarks_', clust,'.RDS'))
            ser_list[[clust]] <- (subser)

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
            return.thresh = Inf, assay = 'RNA', verbose = FALSE)
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
