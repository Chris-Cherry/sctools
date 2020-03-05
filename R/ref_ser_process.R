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
    if(!is.null(ref_ser@meta.data$pct_mito)){
        pct_mt = ref_ser@meta.data$pct_mito
        ser$pct_mito = pct_mt
        scale_vars = c(scale_vars, 'pct_mito')
    }
    if(!is.null(ref_ser@meta.data$G2M.Score) & !is.null(ref_ser@meta.data$S.Score)){
        G2M.Score = ref_ser@meta.data$G2M.Score
        S.Score = ref_ser@meta.data$S.Score
        scale_vars = c(scale_vars, 'G2M.Score', 'S.Score')
    }
    if(!is.null(ref_ser@meta.data$nCount_RNA)){
        scale_vars = c(scale_vars, 'nCount_RNA')
    }

    ############################################################################
    # IMPLEMENT METHOD FOR LIST OF OTHER GENE SETS
    ############################################################################

    ser = ScaleData(ser, vars.to.regress = scale_vars, verbose = FALSE, features = rownames(ser))
    ser = FindVariableFeatures(ser, dispersion.cutoff = c(1, Inf), verbose = FALSE)
    ser = FindNeighbors(ser, reduction = 'mnn', verbose = FALSE)
    ser = FindClusters(ser, resolution = res, verbose = FALSE)
    ser = RunUMAP(ser, reduction = 'mnn', dims = 1:50, verbose = FALSE)
    phate = phateR::phate(ser@reductions$mnn@cell.embeddings, seed = 42, n.jobs = -1, 
        verbose = FALSE)

    ser[['phate']] = CreateDimReducObject(phate$embedding, key = 'PHATE_',
        assay = DefaultAssay(ser))

    # Generate 3d phate    
    phate3d = phateR::phate(ser@reductions$mnn@cell.embeddings, ndim = 3, seed = 42, n.jobs = -1, 
        verbose = FALSE)
    ser[['phate3d']] = CreateDimReducObject(phate3d$embedding, key = 'PHATE3D_',
        assay = DefaultAssay(ser))

    return(ser)
}

