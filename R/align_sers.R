#' Aligns a list of ser objects
#' 
#' This function integrates separate batches of cells typically processed by
#' process_counts_hash.
#' @param sers      A list of Seurat objects to align. The objects should be normalized with variable genes ID'd.
#' @param meta_file Path to a metadata file containing all metadata by sample. Cell barcodes should be sample_name-XXXXXXXX and the metadata file should have sample_name as the first column.
#' @param ref       NULL or integer specifying which ser in the sers list to use as a reference. If NULL then all combinations will be calculated and the best integration pass selected.
#' @param origin    If only 1 seurat object is presented, user have to define parameter to subset the seurat object
#' 
#' @return Integrated Seurat object with mnn components as pca reduction.
#' @import Seurat
#' @import batchelor
#' @import scran
#' @import phateR
#' @import rgl
#' @import SingleCellExperiment
#' @import Matrix
#' @import scater
#' 
#' @export
 

align_sers = function(sers, meta_file = 'metadata.csv', ref = NULL, origin = 'Sample'){
    
    if (length(sers) == 1){
        sers = UpdateSeuratObject(sers)
        sers = SplitObject(sers, split.by = origin)
    }

    genes = c()
    for(ser in sers){
        genes = c(genes, rownames(ser@assays$RNA@counts))
    }
    genes = unique(genes)

    min_cells = Inf
    for(ser in sers){
        if(ncol(ser) < min_cells){min_cells = ncol(ser)}
    }

    sers = lapply(sers, function(ser){
        counts = ser@assays$RNA@counts
        new_counts = matrix(0, ncol = ncol(ser), nrow = length(genes))
        rownames(new_counts) = genes
        colnames(new_counts) = colnames(ser)
        new_counts[rownames(counts),] = matrix(counts)
        ser = Seurat::CreateSeuratObject(new_counts, meta.data = ser[[]])
        ser = Seurat::NormalizeData(ser, verbose = FALSE)
        ser = Seurat::FindVariableFeatures(ser, verbose = FALSE)
    })

    # Convert Seurat objects to SingleCellExperiment objects
    sce = list()
    for(i in 1:length(sers)){
        sce[i] <- list(as.SingleCellExperiment(sers[[i]]))
        sce[[i]] <- computeSumFactors(sce[[i]])
        sce[[i]] <- normalize(sce[[i]])
    }
        
    # Initate a single sce object of the combined data
    for(i in 1:(length(sce) - 1)){
        if (i == 1){
            sce_counts <- cbind(sce[[i]]@assays@data@listData$counts, sce[[i+1]]@assays@data@listData$counts)
            sce_logcount <- cbind(sce[[i]]@assays@data@listData$counts, sce[[i+1]]@assays@data@listData$logcounts)
            sce_col <- rbind(colData(sce[[i]]), colData(sce[[i+1]]))
        }
        else {
            sce_counts <- cbind(sce_counts, sce[[i+1]]@assays@data@listData$counts)
            sce_logcount <- cbind(sce_logcount, sce[[i+1]]@assays@data@listData$logcounts)
            sce_col <- rbind(sce_col, colData(sce[[i+1]]))
        }
    }

    sce_merge <- SingleCellExperiment(assays = list(counts = sce_counts, logcounts = sce_logcount), rowData = rowData(sce[[1]]),
        colData = sce_col)

    # Multiple batch correction
    rescale <- multiBatchNorm(sce_merge, batch = eval(parse(text = paste0("sce_merge@colData@listData$", origin))))

    mnn_out <- fastMNN(rescale, batch = eval(parse(text = paste0("sce_merge@colData@listData$", origin))), 
    subset.row = rownames(sce[[1]]), k = 20, d = 50, BNPARAM = BiocNeighbors::AnnoyParam())

    reducedDim(sce_merge, "mnn") <- reducedDim(mnn_out, "corrected")

    ser = as.Seurat(sce_merge)
    
    metadata = read.table(meta_file, sep = ',', header = TRUE, row.names = 1)
    ser_samples = sapply(colnames(ser), function(x){
        strsplit(x, '-', fixed = TRUE)[[1]][1]
    })
    ser_metas = metadata[ser_samples,]
    rownames(ser_metas) = names(ser_samples)
    ser = Seurat::AddMetaData(ser, ser_metas)

    return(ser)
}
