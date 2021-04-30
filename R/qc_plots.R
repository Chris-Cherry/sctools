#' Probs need documentation here; or should mush into the process_counts_hash function

qc_plots <- function(ser, dir = "0.1_qc", umi_thresh = 500, feature_thresh = 0.001,
    mt_handle = "^mt-"){
    
    ifelse(!dir.exists(file.path(dir)), dir.create(file.path(dir)),FALSE)
    setwd(dir)

    mt_genes = grep(mt_handle, rownames(ser@assays$RNA@counts))
    dat = Seurat::GetAssayData(ser, slot = 'counts', assay = 'RNA')
    pct_mt = Matrix::colSums(dat[mt_genes,])/ser$nCount_RNA
    ser$pct_mt = pct_mt

    pdf('histograms.pdf')   
        hist(ser$nCount_RNA, breaks = 200, 
            main = paste0('Total UMI count for all cells'))
        lines(x = c(umi_thresh, umi_thresh), y = c(0, 1000), col = 'red')
        hist(ser$nFeature_RNA, breaks = 200,
            main = paste0('Total feature count for all cells'))
        lines(x = c(feature_thresh*ncol(ser), feature_thresh*ncol(ser)), y = c(0, 1000), col = 'red')
    dev.off()

    # Create scatter plots for cells feature and umi counts with thresholds.
    pdf('scatterplot.pdf')
        plot(ser$nCount_RNA, ser$nFeature_RNA, cex = .1, 
            main = 'UMI vs Feature count for all cells')
        lines(x = c(umi_thresh, umi_thresh), y = c(-10000, 10000), col = 'red')
        lines(x = c(-100000, 100000), y = c(feature_thresh*ncol(ser), feature_thresh*ncol(ser)), col = 'red')

        plot(ser$nCount_RNA, ser$nFeature_RNA, cex = .1, 
            main = 'UMI vs Feature count for all cells', xlim = c(0, 10000), 
            ylim = c(0, 2000))
        lines(x = c(umi_thresh, umi_thresh), y = c(-10000, 10000), col = 'red')
        lines(x = c(-100000, 100000), y = c(feature_thresh*ncol(ser), feature_thresh*ncol(ser)), col = 'red')
    dev.off()

    pdf('vln_plots.pdf')
        print(Seurat::VlnPlot(ser, c('nCount_RNA', 'nFeature_RNA', 'pct_mt'), pt.size = 0, ncol = 3))
    dev.off()
}