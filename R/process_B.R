#' This function integrates B cell data vdj sequence into the seurat object
#' 
#' @param ser         Seurat object to process
#' @param B_dir         Optional directory to B cell BCR
#' 
#' @import Seurat
#' @importFrom utils read.table
#' @importFrom data.table as.data.table
#' @importFrom data.table data.table
#' @return Outputs a Seurat with B cell integration
#' 

process_B <- function(ser, B_dir){

    Bcell = read.table(B_dir, sep = ',', stringsAsFactors = FALSE)$V1
    Bcell = Bcell[2:length(Bcell)]
    Bcell = sapply(Bcell, function(x){strsplit(x, '-', fixed = TRUE)[[1]][1]})
    names(Bcell) = Bcell
    Bcell_chain = read.table(B_dir, sep = ',', stringsAsFactors = FALSE)$V6
    Bcell_chain = Bcell_chain[2:length(Bcell_chain)]
    Bcell_seq = read.table(B_dir, sep = ',', stringsAsFactors = FALSE)$V13
    Bcell_seq = Bcell_seq[2:length(Bcell_seq)]
    Bcell_dt = data.table(Bcell, Bcell_chain, Bcell_seq)
    Bcell_mx = as.data.table(matrix('NA', ncol = 7, nrow = nrow(unique(Bcell_dt[,1]))), stringAsFactors = FALSE)
    colnames(Bcell_mx) = c("Bcell", "IGH", "IGK", "TGL", "IGH_seq", "IGK_seq", "IGL_seq")
    rownames(Bcell_mx) = unlist(unique(Bcell_dt[,1]))
    j = 1
    dups = ""

    # Organize the T cell TCR matrix
    for (i in 1:length(rownames(Bcell_dt))){
        if (is.na(match(Bcell_dt[i, 1], dups))){
            dups = Bcell_dt[i, 1]
            Bcell_mx[j,1] = Bcell_dt[i,1]
            if (Bcell_dt[i,2] == "IGH"){
                Bcell_mx[j,2] = Bcell_dt[i,2]
                Bcell_mx[j,5] = Bcell_dt[i,3]
            }  
            else if(Bcell_dt[i,2] == "IGK"){
                Bcell_mx[j,3] = Bcell_dt[i,2]
                Bcell_mx[j,6] = Bcell_dt[i,3]
            }
            else if(Bcell_dt[i,2] == "IGL"){
                Bcell_mx[j,4] = Bcell_dt[i,2]
                Bcell_mx[j,7] = Bcell_dt[i,3]
            }
            j = j + 1

        }
        else{
            if (Bcell_dt[i, 2] == "IGH" & (Bcell_mx[j - 1, 3] == "NA" || Bcell_mx[j - 1, 3] == "None")){
                Bcell_mx[j - 1,2] = Bcell_dt[i,2]
                Bcell_mx[j - 1,5] = Bcell_dt[i,3]
            }  
            else if(Bcell_dt[i,2] == "IGK" & (Bcell_mx[j - 1, 5] == "NA" || Bcell_mx[j - 1, 5] == "None")){
                Bcell_mx[j - 1,3] = Bcell_dt[i,2]
                Bcell_mx[j - 1,6] = Bcell_dt[i,3]
            }
            else if(Bcell_dt[i,2] == "IGL" & (Bcell_mx[j - 1, 5] == "NA" || Bcell_mx[j - 1, 5] == "None")){
                Bcell_mx[j - 1,4] = Bcell_dt[i,2]
                Bcell_mx[j - 1,7] = Bcell_dt[i,3]
            }
        }
    }
    
    # Find matching B cell with TCR seq in the seurat obj
    col = colnames(ser@assays$RNA@counts)
    col = sapply(col, function(x){strsplit(x, '-', fixed = TRUE)[[1]][2]})
    B_ref = as.data.table(matrix('NA', ncol = 7, nrow = length(names(col))))
    for (i in 1:length(col)){
        if (length(which(Bcell_mx[,1] == col[[i]]))){
            B_ref[i,] = Bcell_mx[i,]
        }
        else{
            B_ref[i,1] = col[[i]]
        }
    }
    colnames(B_ref) = colnames(Bcell_mx)
    rownames(B_ref) = colnames(ser@assays$RNA@counts)
    ser = Seurat::AddMetaData(ser, B_ref)

    return(ser)
}