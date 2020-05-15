#' Process B cell clonality
#'
#' This script reads in a counts file from either the DropSeq or 10x 
#' pipeline, converts the genes to a given naming convention (MGI or HGNC) and return a Seurat object 
#'
#' @param ser           Seurat object containing T cell information
#' 
#' @import Seurat
#' @importFrom data.table as.data.table
#' @importFrom plyr compact
#' @import ggplot2
#' 
#' @return Feature plots highlighting clones
#'
#' @export

Bclone_processing <- function(ser){
    
    Bcell_mx = as.data.table(matrix('NA', ncol = 7, nrow = ncol(ser)))
    colnames(Bcell_mx) = c("Bcell", "IGH", "IGK", "IGL", "IGH_seq", "IGK_seq", "IGL_seq")

    Bcell_mx$Bcell = ser@meta.data$Bcell
    Bcell_mx$IGH = ser@meta.data$IGH
    Bcell_mx$IGK = ser@meta.data$IGK
    Bcell_mx$IGL = ser@meta.data$IGL
    Bcell_mx$IGH_seq = ser@meta.data$IGH_seq
    Bcell_mx$IGK_seq = ser@meta.data$IGK_seq
    Bcell_mx$IGL_seq = ser@meta.data$IGL_seq
    rownames(Bcell_mx) = colnames(ser)

    # Find T cell clone
    IGH_IGK_IGL_clone = list()
    IGH_IGK_clone = list()
    IGH_IGL_clone = list()
    single_IGH_clone = list()
    single_IGK_clone = list()
    single_IGL_clone = list()
    j = 1
    k = 1
    m = 1
    n = 1
    for (i in 1:length(rownames(Bcell_mx))){
        if(length(which(Bcell_mx[[i,5]] == Bcell_mx[,5])) > 1 & 
            length(which(Bcell_mx[[i,6]] == Bcell_mx[,6])) > 1 &
            Bcell_mx[[i,5]] != "NA" & Bcell_mx[[i,5]] != "None" &
            Bcell_mx[[i,6]] != "NA" & Bcell_mx[[i,6]] != "None"){
            IGH_IGK_clone[j] <- list(Bcell_mx[which(Bcell_mx[[i,5]] == Bcell_mx[,5]),])
            j = j + 1
        }
        else if(length(which(Bcell_mx[[i,5]] == Bcell_mx[,5])) > 1 & 
            length(which(Bcell_mx[[i,7]] == Bcell_mx[,7])) > 1 &
            Bcell_mx[[i,5]] != "NA" & Bcell_mx[[i,5]] != "None" &
            Bcell_mx[[i,7]] != "NA" & Bcell_mx[[i,7]] != "None"){
            IGH_IGL_clone[j] <- list(Bcell_mx[which(Bcell_mx[[i,5]] == Bcell_mx[,5]),])
            j = j + 1
        }
        else if(length(which(Bcell_mx[[i,5]] == Bcell_mx[,5])) > 1 & Bcell_mx[[i,5]] != "NA" & Bcell_mx[[i,5]] != "None"){
            single_IGH_clone[k] <- list(Bcell_mx[i,])
            k = k + 1
        }
        else if(length(which(Bcell_mx[[i,6]] == Bcell_mx[,6])) > 1 & Bcell_mx[[i,6]] != "NA" & Bcell_mx[[i,6]] != "None"){
            single_IGK_clone[m] <- list(Bcell_mx[i,])
            m = m + 1
        }
        else if(length(which(Bcell_mx[[i,7]] == Bcell_mx[,7])) > 1 & Bcell_mx[[i,7]] != "NA" & Bcell_mx[[i,7]] != "None"){
            single_IGL_clone[n] <- list(Bcell_mx[i,])
            n = n + 1
        }
    }

    single_IGH_clone = compact(unique(single_IGH_clone))
    single_IGH_clone = as.data.frame(do.call(rbind, single_IGH_clone))

    single_IGK_clone = compact(unique(single_IGK_clone))
    single_IGK_clone = as.data.frame(do.call(rbind, single_IGK_clone))

    single_IGL_clone = compact(unique(single_IGL_clone))
    single_IGL_clone = as.data.frame(do.call(rbind, single_IGL_clone))
  
    IGH_IGK_clone = compact(unique(IGH_IGK_clone))
    IGH_IGK_clone = as.data.frame(do.call(rbind, IGH_IGK_clone))

    IGH_IGL_clone = compact(unique(IGH_IGL_clone))
    IGH_IGL_clone = as.data.frame(do.call(rbind, IGH_IGL_clone))

    if (nrow(single_IGH_clone) != 0){
        single_IGH_cell = c()
        for (i in 1:nrow(single_IGH_clone)){
        single_IGH_cell[i] = colnames(ser)[which(single_IGH_clone$Bcell[i] == ser@meta.data$Bcell)]
        }
        p1 <- DimPlot(ser, pt.size = 1, cells.highlight = single_IGH_cell, sizes.highlight = 0.5) +
        scale_color_manual(labels = c("Other cells","IgH clone"), values = c("grey", "darkred")) + theme(legend.text = element_text(size = 10))
    }else{p1 <- c()}
   

    if (nrow(single_IGK_clone) != 0){
        single_IGK_cell = c()
        for (i in 1:nrow(single_IGK_clone)){
            single_IGK_cell[i] = colnames(ser)[which(single_IGK_clone$Bcell[i] == ser@meta.data$Bcell)]
        }
        p2 <- DimPlot(ser, pt.size = 1, cells.highlight = single_IGK_cell , sizes.highlight = 0.5) +
        scale_color_manual(labels = c("Other cells","IgK clone"), values = c("grey", "darkblue")) + theme(legend.text = element_text(size = 10))
    }else {p2 <- c()}
   
    if (nrow(single_IGL_clone) != 0){
        single_IGL_cell = c()
        for (i in 1:nrow(single_IGL_clone)){
            single_IGL_cell[i] = colnames(ser)[which(single_IGL_clone$Bcell[i] == ser@meta.data$Bcell)]
        }
        p3 <- DimPlot(ser, pt.size = 1, cells.highlight = single_IGL_cell, sizes.highlight = 0.5) +
        scale_color_manual(labels = c("Other cells","IgL clone"), values = c("grey", "darkgreen")) + theme(legend.text = element_text(size = 10))
    }else {p3 <- c()}

    if (nrow(IGH_IGK_clone) != 0){
        IGH_IGK_cell = c()
        for (i in 1:nrow(IGH_IGK_clone)){
            IGH_IGK_cell[i] = colnames(ser)[which(IGH_IGK_clone$Bcell[i] == ser@meta.data$Bcell)]
        }
        p4 <- DimPlot(ser, pt.size = 1, cells.highlight = IGH_IGK_cell, sizes.highlight = 0.5) +
        scale_color_manual(labels = c("Other cells","IgH_IgK clone"), values = c("grey", "blue")) + theme(legend.text = element_text(size = 10))
    }else {p4 <- c()}
   
    if (nrow(IGH_IGL_clone) != 0){
        IGH_IGL_cell = c()
        for (i in 1:nrow(IGH_IGL_clone)){
            IGH_IGL_cell[i] = colnames(ser)[which(IGH_IGL_clone$Bcell[i] == ser@meta.data$Bcell)]
        }
        p5 <- DimPlot(ser, pt.size = 1, cells.highlight = IGH_IGL_cell, sizes.highlight = 0.5) +
        scale_color_manual(labels = c("Other cells","IgH_IgL clone"), values = c("grey", "green")) + theme(legend.text = element_text(size = 10))
    }else {p5 <- c()}

    CombinePlots(plots = list(p1,p2,p3,p4,p5), ncol = 2)

}