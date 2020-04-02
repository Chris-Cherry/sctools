#' This file process B cell clonality. 
#' 
#' Parameters:
#' @param ser           Seurat object containing T cell information
#' 
#' This script reads in a counts file from either the DropSeq or 10x 
#' pipeline, converts the genes to a given naming convention (MGI or HGNC) and return a Seurat object 
#' 
#' @import Seurat
#' @import Matrix
#' @import methods
#' @import plyr
#' @import utils
#' @import ggplot2
#' 
#' @export

Bclone_processing <- function(ser){
    
    Bcell_mx = as.data.table(matrix('NA', ncol = 5, nrow = ncol(ser)))
    colnames(Bcell_mx) = c("Bcell", "a_Bchain", "a_Bseq", "b_Bchain", "b_Bseq")

    Bcell_mx$Bcell = ser@meta.data$Bcell
    Bcell_mx$B_chain = ser@meta.data$B_chain
    Bcell_mx$B_seq = ser@meta.data$B_seq
    rownames(Bcell_mx) = colnames(ser)

    # Find T cell clone
    doubleB_clone = list()
    singleB_clone = list()
    j = 1
    k = 1
    for (i in 1:length(rownames(Bcell_mx))){
        if(length(which(Bcell_mx[[i,3]] == Bcell_mx[,3])) > 1 & Bcell_mx[[i,3]] != "NA" & Bcell_mx[[i,3]] != "None"){
            doubleB_clone[j] <- list(Bcell_mx[which(Bcell_mx[[i,3]] == Bcell_mx[,3]),])
            j = j + 1
        }
        else{
            singleB_clone[k] <- list(Bcell_mx[i,])
            k = k + 1
        }
    }

    singleB_clone = compact(unique(singleB_clone))
    singleB_clone = as.data.frame(do.call(rbind, singleB_clone))
    rownames(singleB_clone) = singleB_clone[,1]

    doubleB_clone = compact(unique(doubleB_clone))
    doubleB_clone = as.data.frame(do.call(rbind, doubleB_clone))
    rownames(doubleB_clone) = doubleB_clone[,1]
    
    singleB_cell = c()
    for (i in 1:nrow(singleB_clone)){
        singleB_cell[i] = colnames(ser)[which(singleB_clone$Bcell[i] == ser@meta.data$Bcell)]
    }

    doubleB_cell = c()
    for (i in 1:nrow(doubleB_clone)){
        doubleB_cell[i] = colnames(ser)[which(doubleB_clone$Bcell[i] == ser@meta.data$Bcell)]
    }

    DimPlot(ser, pt.size = 2, cells.highlight = list(singleB_cell, doubleB_cell)) +
    scale_color_manual(labels = c("", "single_clone","double_clone"), values = c("grey", "darkred", "darkblue"))
  
    DimPlot(ser, pt.size = 2, cells.highlight = singleB_cell) +
    scale_color_manual(labels = c("","single_clone"), values = c("grey", "darkred"))
    
    DimPlot(ser, pt.size = 2, cells.highlight = bB_cell) +
    scale_color_manual(labels = c("","double_clone"), values = c("grey", "darkblue"))

}