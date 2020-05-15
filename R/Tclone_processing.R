#' This file process T cell clonality. 
#'
#' This script reads in a counts file from either the DropSeq or 10x 
#' pipeline, converts the genes to a given naming convention (MGI or HGNC) and return a Seurat object 
#' 
#' @param ser           Seurat object containing T cell information
#'
#' @import Seurat
#' @importFrom plyr compact
#' @importFrom data.table as.data.table
#' @return Feature plots highlighting clones
#' @export

Tclone_processing <- function(ser){
    
    Tcell_mx = as.data.table(matrix('NA', ncol = 5, nrow = ncol(ser)))
    colnames(Tcell_mx) = c("Tcell", "a_Tchain", "a_Tseq", "b_Tchain", "b_Tseq")

    Tcell_mx$Tcell = ser@meta.data$Tcell
    Tcell_mx$a_Tchain = ser@meta.data$a_Tchain
    Tcell_mx$a_Tseq = ser@meta.data$a_Tseq
    Tcell_mx$b_Tchain = ser@meta.data$b_Tchain
    Tcell_mx$b_Tseq = ser@meta.data$b_Tseq
    rownames(Tcell_mx) = colnames(ser)

    # Find T cell clone
    aT_clone = list()
    bT_clone = list()
    doubleT_clone = list()
    singleT_clone = list()
    j = 1
    k = 1
    m = 1
    n = 1
    for (i in 1:length(rownames(Tcell_mx))){
        if(length(which(Tcell_mx[[i,3]] == Tcell_mx[,3])) > 1 & length(which(Tcell_mx[[i,5]] == Tcell_mx[,5])) > 1 &
            Tcell_mx[[i,3]] != "NA" & Tcell_mx[[i,3]] != "None" & Tcell_mx[[i,5]] != "NA" & Tcell_mx[[i,5]] != "None"
            & length(intersect(which(Tcell_mx[[i,3]] == Tcell_mx[,3]), which(Tcell_mx[[i,5]] == Tcell_mx[,5]))) > 1){
            doubleT_clone[j] <- list(Tcell_mx[intersect(which(Tcell_mx[[i,3]] == Tcell_mx[,3]), which(Tcell_mx[[i,5]] == Tcell_mx[,5])),])
            j = j + 1
        }
        else if(length(which(Tcell_mx[[i,3]] == Tcell_mx[,3])) > 1 & Tcell_mx[i,3] != "NA" & Tcell_mx[i,3] != "None"){
            aT_clone[k] <- list(Tcell_mx[which(Tcell_mx[[i,3]] == Tcell_mx[,3]),])
            k = k + 1
        }
        else if(length(which(Tcell_mx[[i,5]] == Tcell_mx[,5])) > 1 & Tcell_mx[i,5] != "NA" & Tcell_mx[i,5] != "None"){
            bT_clone[m] <- list(Tcell_mx[which(Tcell_mx[[i,5]] == Tcell_mx[,5]),])
            m = m + 1
        }
        else{
            singleT_clone[n] <- list(Tcell_mx[i,])
            n = n + 1
        }
    }

    aT_clone = compact(unique(aT_clone))
    if (nrow(aT_clone[[1]]) == 0){
        aT_clone = aT_clone[2:length(aT_clone)]
    }
    aT_clone = as.data.frame(do.call(rbind, aT_clone))
    rownames(aT_clone) = aT_clone[,1]

    bT_clone = compact(unique(bT_clone))
    if (nrow(bT_clone[[1]]) == 0){
        bT_clone = bT_clone[2:length(bT_clone)]
    }
    bT_clone = as.data.frame(do.call(rbind, bT_clone))
    rownames(bT_clone) = bT_clone[,1]

    singleT_clone = compact(unique(singleT_clone))
    singleT_clone = as.data.frame(do.call(rbind, singleT_clone))
    rownames(singleT_clone) = singleT_clone[,1]

    doubleT_clone = compact(unique(doubleT_clone))
    doubleT_clone = as.data.frame(do.call(rbind, doubleT_clone))
    rownames(doubleT_clone) = doubleT_clone[,1]

    if (nrow(aT_clone) !=0){
        aT_cell = c()
        for (i in 1:nrow(aT_clone)){
            aT_cell[i] = colnames(ser)[which(aT_clone$Tcell[i] == ser@meta.data$Tcell)]
            p1 <- DimPlot(ser, pt.size = 1, cells.highlight = aT_cell, sizes.highlight = 0.5) +
            scale_color_manual(labels = c("Other cells","a_clone"), values = c("grey", "darkred")) + theme(legend.text = element_text(size = 10))
        }
    }else{p1<-c()}
    
    if (nrow(bT_clone) !=0){
        bT_cell = c()
        for (i in 1:nrow(bT_clone)){
            bT_cell[i] = colnames(ser)[which(bT_clone$Tcell[i] == ser@meta.data$Tcell)]
        }
        p2 <- DimPlot(ser, pt.size = 1, cells.highlight = bT_cell, sizes.highlight = 0.5) +
        scale_color_manual(labels = c("Other cells","b_clone"), values = c("grey", "darkblue")) + theme(legend.text = element_text(size = 10))
    }else{p2<-c()}


    if (nrow(doubleT_clone) !=0){
        doubleT_cell = c()
        for (i in 1:nrow(doubleT_clone)){
            doubleT_cell[i] = colnames(ser)[which(doubleT_clone$Tcell[i] == ser@meta.data$Tcell)]
        }
        p3 <- DimPlot(ser, pt.size = 1, cells.highlight = doubleT_cell, sizes.highlight = 0.5) +
        scale_color_manual(labels = c("Other cells","ab_clone"), values = c("grey", "darkgreen")) + theme(legend.text = element_text(size = 10))
    }else{p3<-c()}
     
    CombinePlots(plots = list(p1, p2, p3))

}