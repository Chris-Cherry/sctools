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

    if(min_cells < 200){
        k.filter = min_cells*.5
    }else{k.filter = 200}
    if(min_cells < 30){
        k.score = min_cells
    }else{k.score = 30}
    if(min_cells < 5){
        k.anchor = min_cells
    }else{k.anchor = 5}

    features = Seurat::SelectIntegrationFeatures(sers, verbose = FALSE)
    anchors = Seurat::FindIntegrationAnchors(sers, verbose = FALSE, 
        k.filter = k.filter, k.score = k.score, k.anchor = k.anchor)
    ser = Seurat::IntegrateData(anchors, features.to.integrate = genes, verbose = FALSE)
    Seurat::DefaultAssay(ser) = 'integrated'
    
    metadata = read.table(meta_file, sep = ',', header = TRUE, row.names = 1)
    ser_samples = paste0(sapply(colnames(ser), function(x){
        strsplit(x, '_', fixed = TRUE)[[1]][2]
    }), '_', sapply(colnames(ser), function(x){
        strsplit(x, '_', fixed = TRUE)[[1]][3]
    }))
    ser_metas = metadata[ser_samples,]
    rownames(ser_metas) = colnames(ser)
    ser = AddMetaData(ser, ser_metas)

    return(ser)
}
