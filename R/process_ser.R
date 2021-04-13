#' Processes counts into a Seurat object prepared for alignment.
#' 
#' Reads in a blank ser object (usually from 2) and processes
#' with a traditional Seurat pipeline. By default will scale both RNA and
#' the default assay if not RNA but will perform PCA on default assay.
#'
#' @param ser           Seurat object to process.
#' @param mt_handle     Regex used to identify mitochondrial genes for scaling. If left blank mt gene percent will not be used to scale.
#' @param mt_thresh     Max percent of umis from mitochondrial genes for cells to be used. If mt_handle isn't used this won't have any functionality.
#' @param feat_thresh   Proportion of cells a gene must be expressed in to be included
#' @param scale_umi     Whether or not to scale on total UMI count.
#' @param g2m_genes     Genes to use for g2m scoring and scaling. If left blank cell cycle scoring and scaling will not be done.
#' @param s_genes       Genes to use for s scoring and scaling. If left blank cell cycle scoring and scaling will not be done.
#' @param res           Resolution for clustering.
#' @param other_sets    A named list of gene sets to be used similar to percent mt for scoring and scaling. Names will appear in metadata.
#' @param ref_ser       A processed reference Seurat object used to as reference for cell selection.
#' @param scale_vars    Other features to use for scaling. 
#' @param phate         Boolean for whether to create phate embeddings
#' @param verbose       Boolean for whether to run verbose versions of functions
#'
#' @import Seurat
#' @importFrom phateR phate
#' @importFrom Matrix rowSums
#' @importFrom Matrix colSums
#' @return Outputs a processed Seurat outputs (UMAP, Phate) 
#' @export

process_ser <- function(ser, mt_handle = NULL, mt_thresh = .1, feat_thresh = 0.001, scale_umi = TRUE, 
    g2m_genes = NULL, s_genes = NULL, res = .8, other_sets = NULL, ref_ser = NULL,
    scale_vars = NULL, phate = FALSE, verbose = FALSE){

    ser = Seurat::UpdateSeuratObject(ser)

    feat_sums = Matrix::rowSums(Seurat::GetAssayData(ser, slot = 'counts', assay = 'RNA') != 0)
    feat_keep = names(feat_sums)[which(feat_sums > ncol(ser)*feat_thresh)]
    ser = subset(ser, features = feat_keep)
    
    if(!is.null(ref_ser)){
        ser = ser[ , intersect(colnames(ser), colnames(ref_ser))]
    }

    if(is.null(scale_vars)){scale_vars = c()}

    if(!is.null(mt_handle)){
        mt_genes = grep(mt_handle, rownames(ser@assays$RNA@counts))
        dat = Seurat::GetAssayData(ser, slot = 'counts', assay = 'RNA')
        pct_mt = Matrix::colSums(dat[mt_genes,])/ser$nCount_RNA
        ser$pct_mt = pct_mt
        keep = names(pct_mt)[which(pct_mt <= mt_thresh)]
        ser = subset(ser, cells = keep)
        scale_vars = c(scale_vars, 'pct_mt')
    }

    if(!is.null(s_genes) & !is.null(g2m_genes)){
        ser = Seurat::CellCycleScoring(ser, readLines(s_genes), readLines(g2m_genes), assay = "RNA")
        scale_vars = c(scale_vars, 'G2M.Score', 'S.Score')
    }

    if(scale_umi){
        scale_vars = c(scale_vars, 'nCount_RNA')
    }
    
    if(Seurat::DefaultAssay(ser) != 'RNA'){
        ser = Seurat::ScaleData(ser, vars.to.regress = scale_vars, assay = 'RNA', 
            verbose = verbose, features = rownames(ser@assays$RNA@data))
    }
    ser = Seurat::ScaleData(ser, vars.to.regress = scale_vars, verbose = verbose, features = rownames(ser))
    ser = Seurat::FindVariableFeatures(ser, verbose = verbose)
    ser = Seurat::RunPCA(ser, npcs =50,verbose = verbose)
    ser = Seurat::FindNeighbors(ser, reduction = "pca", verbose = verbose)
    ser = Seurat::FindClusters(ser, resolution = res, verbose = verbose)
    ser = Seurat::RunUMAP(ser, reduction = "pca", dims = 1:50, verbose = verbose)
    
    if (phate){
        phate = phateR::phate(ser@reductions$pca@cell.embeddings, seed = 42, n.jobs = -1, 
            verbose = verbose)

        ser[['phate']] = Seurat::CreateDimReducObject(100*phate$embedding, key = 'PHATE_',
            assay = DefaultAssay(ser))

        # Generate 3d phate    
        phate3d = phateR::phate(ser@reductions$pca@cell.embeddings, ndim = 3, seed = 42, n.jobs = -1, 
            verbose = verbose)
        ser[['phate3d']] = Seurat::CreateDimReducObject(phate3d$embedding, key = 'PHATE3D_',
            assay = DefaultAssay(ser))
    }

    return(ser)
}
