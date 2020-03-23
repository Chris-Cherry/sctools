#' Processes counts into a Seurat object prepared for alignment. 
#' This file also seperated.Adds hashing information to a Seurat object and adds hash id to barcodes.
#' Parameters:
#' @param directory     directory containing the output files of align_dropseq_2.sh 
#'                      or cell ranger
#' @param from_gene     "ENSG" or "ENSMUSG" or "HGNC" or "MGI" 
#' @param to_gene       "MGI" or "HGNC"
#' @param sample_id     Optional experimental ID to prefix cell names
#' @param hash_dir      Directory containing hash counts files
#' @param hash_meta     csv file containing conversion of hash ID to sample name
#' @param type          "DropSeq" or "10x"
#' @param T_dir         Optional directory to T cell TCR
#' @param B_dir         optional directory to B cell BCR
#' @param umi_thresh    Minimum total umi count for cells
#' 
#' This script reads in a counts file from either the DropSeq or 10x 
#' pipeline, converts the genes to a given naming convention (MGI or HGNC) and return a Seurat object 
#' 
#' @import Seurat
#' @import Matrix
#' @import methods
#' @import plyr
#' @import utils
#' @import data.table
#' @return Outputs a Seurat file prepared to rPCA anchoring
#' @export

process_counts_hash <- function(directory, from_gene, to_gene, hash_dir = NULL, 
    hash_meta = NULL, sample_id = NULL, type = 'DropSeq', T_dir = NULL, 
    B_dir = NULL, umi_thresh = 500){



    if(type == 'DropSeq'){
        prefix = tail(strsplit(directory, '/', fixed = TRUE)[[1]], n = 1)
        file = paste0(directory, '/', prefix, '_DGE.txt.gz')
        counts =  read.table(gzfile(file), header = TRUE, row.names = 1)
    }
    
    if(type == '10x'){
        counts = readMM(gzfile(paste0(directory, 
            '/filtered_feature_bc_matrix/matrix.mtx.gz')))
        rows = read.table(gzfile(paste0(directory,
            '/filtered_feature_bc_matrix/features.tsv.gz')), 
            stringsAsFactors = FALSE)$V1
        rows = sapply(rows, function(x){
            strsplit(x, '-', fixed = TRUE)[[1]][1]
        })
        names(rows) = c()
        rownames(counts) = rows

        cols = read.table(gzfile(paste0(directory,
            '/filtered_feature_bc_matrix/barcodes.tsv.gz')), 
            stringsAsFactors = FALSE)$V1
        cols = sapply(cols, function(x){
            strsplit(x, '-', fixed = TRUE)[[1]][1]
        })
        names(cols) = c()
        colnames(counts) = cols
    }

    if(!is.null(sample_id)){
        colnames(counts) = paste(sample_id, colnames(counts), sep = '-')
    }

    c_sum = colSums(counts)
    c_take = which(c_sum >= umi_thresh)
    counts = counts[,names(c_take)]
    counts = as(as.matrix(counts), "sparseMatrix")

    # Convert genes
    if(from_gene != to_gene){
        conv = convert_genes(rownames(counts), from = from_gene, to = to_gene)
        counts = counts[-which(is.na(match(rownames(counts), conv[,1]))),]
        conv = conv[match(rownames(counts), conv[,1]),]
        rm_id = which(conv[,2] == '')
        if (!identical(rm_id, integer(0))){
            counts = counts[-rm_id,]
            conv = conv[-rm_id,]
        }
        
        dups = unique(conv[which(duplicated(conv[,2])),2])
        for(dup in dups){
            rows = which(conv[,2] == dup)
            counts[rows[1],] = colSums(counts[rows,])
            counts = counts[-rows[-1],]
            conv = conv[-rows[-1],]
        }
        
        rownames(counts) = conv[,2]
    }

    ser = CreateSeuratObject(counts)
    
     if(!is.null(hash_dir) & !is.null(hash_meta)){
        if(type == 'DropSeq'){
        # IMPLEMENT DROPSEQ HASHING PROTOCOL
        } else if(type == '10x'){
            counts = readMM(gzfile(paste0(hash_dir, 
                '/filtered_feature_bc_matrix/matrix.mtx.gz')))
            rows = read.table(gzfile(paste0(hash_dir,
                '/filtered_feature_bc_matrix/features.tsv.gz')), 
                stringsAsFactors = FALSE)$V1
            rows = sapply(rows, function(x){
                strsplit(x, '-', fixed = TRUE)[[1]][1]
            })
            names(rows) = c()
            rownames(counts) = rows

            cols = read.table(gzfile(paste0(hash_dir,
                '/filtered_feature_bc_matrix/barcodes.tsv.gz')), 
                stringsAsFactors = FALSE)$V1
            cols = sapply(cols, function(x){
                strsplit(x, '-', fixed = TRUE)[[1]][1]
            })
            names(cols) = c()
            colnames(counts) = cols
        }

        cells = intersect(colnames(counts), colnames(ser))
        counts = counts[, cells]
        ser = subset(ser, cells = cells)
        ser[['hash']] = CreateAssayObject(counts)

        ser = NormalizeData(ser, assay = 'hash', normalization.method = 'CLR')
        ser = HTODemux(ser, assay = 'hash', positive.quantile = .99)
        singlets = names(ser$hash_classification.global)[
            which(ser$hash_classification.global == 'Singlet')]
        ser = subset(ser, cells = singlets)

        hash_conv = read.table(hash_meta, sep = ',', row.names = 1)
        ser[['Sample']] = hash_conv[ser$hash_maxID,]
        meta_drop = match(c('hash_classification', 'hash_classification.global', 
            'hash.ID', 'hash_margin', 'hash_secondID', 'hash_maxID', 
            'nFeature_hash', 'nCount_hash'), colnames(ser@meta.data))
        ser@meta.data = ser@meta.data[,-meta_drop]

        if(type == '10x'){
            nn = sapply(colnames(ser), function(x){
                strsplit(x, '-', fixed = TRUE)[[1]][1]
            })
            names(nn) = c()
            ser = RenameCells(ser, new.names = nn)
        }

        ser = RenameCells(ser, new.names = paste0(ser$Sample, '-', colnames(ser)))
    }
    
     # Matching T with TCR and BCR sequencing
    if (!is.null(T_dir)){
        ser = process_T(ser, T_dir)
    }

    if (!is.null(B_dir)){
        ser = process_B(ser, B_dir)
    }

    return(ser)
}