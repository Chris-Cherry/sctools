#' Processes counts into a Seurat object prepared for alignment.
#' Parameters:
#' @param ser           Seurat object to process.
#' @param mt_handle     Regex used to identify mitochondrial genes for scaling. If 
#'                      left blank mt gene % will not be used to scale.
#' @param scale_umi     Whether or not to scale on total UMI count
#' @param g2m_genes     Genes to use for g2m scoring and scaling. If left blank
#'                      cell cycle scoring and scaling will not be done.
#' @param s_genes       Genes to use for s scoring and scaling. If left blank
#'                      cell cycle scoring and scaling will not be done.
#' @param res           Resolution for clustering
#' @param other_sets    A named list of gene sets to be used similar to %mt for 
#'                      scoring and scaling. Names will appear in metadata.
#' @param ref_ser       A processed reference Seurat object used to as reference
#'                      for cell selection.
#'
#' 
#' Reads in a blank ser object (usually from align_sers.R) and processes
#' with a traditional Seurat pipeline. By default will scale both RNA and
#' the default assay if not RNA but will perform PCA on default assay.
#' 
#' @import Seurat
#' @return  Outputs a processed Seurat outputs (PCA, UMAP, Phate) 
#' @export

process_ser <- function(ser, mt_handle = NULL, scale_umi = TRUE, 
    g2m_genes = NULL, s_genes = NULL, res = .8, other_sets = NULL, ref_ser = NULL){

    feat_sums = Matrix::rowSums(Seurat::GetAssayData(ser, slot = 'counts', assay = 'RNA') != 0)
    feat_keep = names(feat_sums)[which(feat_sums > ncol(ser)*.001)]
    ser = subset(ser, features = feat_keep)
    
    if(!is.null(ref_ser)){
        ser = ser[ , intersect(colnames(ser), colnames(ref_ser))]
    }
    
    if(Seurat::DefaultAssay(ser) != 'integrated'){
        ser = Seurat::NormalizeData(ser, verbose = FALSE)
    }

    ser = Seurat::FindVariableFeatures(ser, verbose = FALSE)
    
    scale_vars = c()
    if(!is.null(mt_handle)){
        mt_genes = grep(mt_handle, rownames(ser))
        dat = Seurat::GetAssayData(ser, slot = 'counts', assay = 'RNA')
        pct_mt = Matrix::colSums(dat[mt_genes,])/ser$nCount_RNA
        ser$pct_mt = pct_mt
        keep = names(pct_mt)[which(pct_mt <= .1)]
        ser = subset(ser, cells = keep)
        scale_vars = c(scale_vars, 'pct_mt')
    }
    if(!is.null(s_genes) & !is.null(g2m_genes)){
        ser = Seurat::CellCycleScoring(ser, s_genes, g2m_genes)
        scale_vars = c(scale_vars, 'G2M.Score', 'S.Score')
    }
    if(scale_umi){
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

