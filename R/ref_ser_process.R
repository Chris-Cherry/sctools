#' Processes counts into a Seurat object prepared for alignment based on reference object.
#' Parameters:
#' @param ser           Seurat object to process.
#' @param res           Resolution for clustering
#' @param ref_ser       A processed reference Seurat object used to obtain metadata
#'                      for scaling and cell selection.
#'
#' Reads in a blank ser object (usually from align_sers.R) and processes
#' with a traditional Seurat pipeline. All scaling will be performed using
#' metadata from the reference Seurat object
#' 
#' @import  Seurat
#' @import  Matrix
#' @import  phateR
#' @return  Outputs a processed Seurat outputs (PCA, UMAP, Phate) 
#' @export

ref_ser_process <- function(ser, res = .8, ref_ser = NULL){

    feat_sums = Matrix::rowSums(Seurat::GetAssayData(ser, slot = 'counts', assay = 'RNA') != 0)
    feat_keep = names(feat_sums)[which(feat_sums > ncol(ser)*.001)]
    ser = subset(ser, features = feat_keep)

    if(!is.null(ref_ser)){
    # Keep same cells as reference
    ser = ser[ , intersect(colnames(ser), colnames(ref_ser))]  
    }

    if(Seurat::DefaultAssay(ser) != 'integrated'){
        ser = Seurat::NormalizeData(ser, verbose = FALSE)
    }

    ser = Seurat::FindVariableFeatures(ser, verbose = FALSE)
    
    scale_vars = c()
    if(exists(ref_ser@meta.data$pct_mt)){
        pct_mt = ref_ser@meta.data$pct_mt
        ser$pct_mt = pct_mt
        scale_vars = c(scale_vars, 'pct_mt')
    }
    if(exists(ref_ser@meta.data$G2M.Score) & exists(ref_ser@meta.data$S.Score)){
        G2M.Score = ref_ser@meta.data$G2M.Score
        S.Score = ref_ser@meta.data$S.Score
        scale_vars = c(scale_vars, 'G2M.Score', 'S.Score')
    }
    if(exists(ref_ser@meta.data$nCount_RNA)){
        scale_vars = c(scale_vars, 'nCount_RNA')
    }

    ############################################################################
    # IMPLEMENT METHOD FOR LIST OF OTHER GENE SETS
    ############################################################################

    if(Seurat::DefaultAssay(ser) != 'RNA'){
        ser = Seurat::ScaleData(ser, vars.to.regress = scale_vars, assay = 'RNA', 
            verbose = FALSE, features = rownames(ser))
    }
    ser = Seurat::ScaleData(ser, vars.to.regress = scale_vars, verbose = FALSE, features = rownames(ser))
    ser = Seurat::RunPCA(ser, npcs = 50, verbose = FALSE)
    ser = Seurat::FindNeighbors(ser, verbose = FALSE)
    ser = Seurat::FindClusters(ser, resolution = res, verbose = FALSE)
    ser = Seurat::RunUMAP(ser, dims = 1:50, verbose = FALSE)
    phate = phateR::phate(ser@reductions$pca@cell.embeddings, seed = 42, n.jobs = -1, 
        verbose = FALSE)
    ser[['phate']] = Seurat::CreateDimReducObject(phate$embedding, key = 'PHATE_',
        assay = Seurat::DefaultAssay(ser))

    # Generate 3d phate    
    phate3d = phateR::phate(ser@reductions$pca@cell.embeddings, ndim = 3, seed = 42, n.jobs = -1, 
        verbose = FALSE)
    ser[['phate3d']] = Seurat::CreateDimReducObject(phate3d$embedding, key = 'PHATE3D_',
        assay = Seurat::DefaultAssay(ser))

    return(ser)
}

