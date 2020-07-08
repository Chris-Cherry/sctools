#' Processes_T: This function integrate T cell data with alpha and beta vdj sequence into the seurat object
#' 
#' @param ser           Seurat object to process
#' @param T_dir           Optional directory to T cell TCR
#' 
#' 
#' @import Seurat
#' @importFrom utils read.table
#' @importFrom data.table as.data.table
#' @return Outputs a Seurat with T cell integration
#' 

process_T <- function(ser, T_dir){

    Tcell = read.table(T_dir, sep = ',', stringsAsFactors = FALSE)$V1
    Tcell = Tcell[2:length(Tcell)]
    Tcell = sapply(Tcell, function(x){strsplit(x, '-', fixed = TRUE)[[1]][1]})
    names(Tcell) = Tcell
    Tcell_chain = read.table(T_dir, sep = ',', stringsAsFactors = FALSE)$V6
    Tcell_chain = Tcell_chain[2:length(Tcell_chain)]
    Tcell_seq = read.table(T_dir, sep = ',', stringsAsFactors = FALSE)$V13
    Tcell_seq = Tcell_seq[2:length(Tcell_seq)]
    Tcell_dt = data.table(Tcell, Tcell_chain, Tcell_seq)
    Tcell_mx = as.data.table(matrix('NA', ncol = 5, nrow = nrow(unique(Tcell_dt[,1]))), stringAsFactors = FALSE)
    colnames(Tcell_mx) = c("Tcell", "a_Tchain", "a_Tseq", "b_Tchain", "b_Tseq")
    rownames(Tcell_mx) = unlist(unique(Tcell_dt[,1]))
    j = 1
    dups = ""

    # Organize the T cell TCR matrix
    for (i in 1:length(rownames(Tcell_dt))){
        if (is.na(match(Tcell_dt[i, 1], dups))){
            dups = Tcell_dt[i, 1]
            Tcell_mx[j,1] = Tcell_dt[i,1]
            if (Tcell_dt[i,2] == "TRA"){
                Tcell_mx[j,2] = Tcell_dt[i,2]
                Tcell_mx[j,3] = Tcell_dt[i,3]
            }  
            else if(Tcell_dt[i,2] == "TRB"){
                Tcell_mx[j,4] = Tcell_dt[i,2]
                Tcell_mx[j,5] = Tcell_dt[i,3]
            }
            j = j + 1
        }
        else{
            if (Tcell_dt[i, 2] == "TRA" & (Tcell_mx[j - 1, 3] == "NA" || Tcell_mx[j - 1, 3] == "None")){
                Tcell_mx[j - 1,2] = Tcell_dt[i,2]
                Tcell_mx[j - 1,3] = Tcell_dt[i,3]
            }  
            else if(Tcell_dt[i,2] == "TRB" & (Tcell_mx[j - 1, 5] == "NA" || Tcell_mx[j - 1, 5] == "None")){
                Tcell_mx[j - 1,4] = Tcell_dt[i,2]
                Tcell_mx[j - 1,5] = Tcell_dt[i,3]
            }
        }
    }
    
    # Find matching T cell with TCR seq in the seurat obj
    col = colnames(ser@assays$RNA@counts)
    col = sapply(col, function(x){strsplit(x, '-', fixed = TRUE)[[1]][2]})
    T_ref = as.data.table(matrix('NA', ncol = 5, nrow = length(names(col))))
    for (i in 1:length(col)){
        if (length(which(Tcell_mx[,1] == col[[i]]))){
            T_ref[i,] = Tcell_mx[which(Tcell_mx[,1] == col[[i]]),]
        }
        else{
            T_ref[i,1] = col[[i]]
        }
    }
    colnames(T_ref) = colnames(Tcell_mx)
    rownames(T_ref) = colnames(ser@assays$RNA@counts)
    ser = AddMetaData(ser, T_ref)
   
    return(ser)
}