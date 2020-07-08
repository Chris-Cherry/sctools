#' Generate flashing and spinning 3d pseudotime plots
#' 
#' @param ser               Seurat object to process
#' @param curve_dir         Input curve directory
#' @param pseudotime_dir    Input pseudotime directory
#' @param bg_color          Background color    
#' @param frame             Number of frames
#' @param out_dir           Output directory
#' @param subset            A subset of clusters that user want to see
#' @import Seurat 
#' @import rgl
#' @importFrom grDevices hcl
#' @return Makes png files that can be assembled into 3d moving plots
#' @export

pseudotime_plot3d <- function(ser, curve_dir, pseudotime_dir, subset = NULL, out_dir = '5_pseudotime_3d', bg_color = 'white', frame = 500){

    # If user want to subset

    if(!is.null(subset)){
        p3 = subset(ser, idents = subset)
        p3 = p3@reductions$phate3d@cell.embeddings
    }
    else{
        p3 = ser@reductions$phate3d@cell.embeddings
    }


    curves = readRDS(curve_dir)
    pseudo = readRDS(pseudotime_dir)

    gg_low_alpha <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, alpha = .1, c = 100)[1:n]
    }
    cols = gg_low_alpha(18)
    names(cols) = 0:17
    dim_cell_cols = cols[Idents(ser)]
    names(dim_cell_cols) = names(Idents(ser))

    gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
    }
    cols = gg_color_hue(18)
    names(cols) = 0:17
    cell_cols = cols[Idents(ser)]
    names(cell_cols) = names(Idents(ser))

    dir.create(out_dir)
    dir.create(paste0(out_dir,'/wflashing'))
    dir.create(paste0(out_dir,'/wspinning'))
    # Flashing and spinning
    for(i in 100:1000){
        pct = i%%100
        alphas = rep(0.25, nrow(p3))
        names(alphas) = rownames(p3)
        for(j in 1:length(pseudo)){
            times = sort(pseudo[[j]])
            times = times[-c(1:400)]
            window_min = max(times)*(pct-5)/100
            window_max = max(times)*(pct+5)/100
            tar_cells = names(times[which(times>window_min & times<window_max)])
            n_cells = length(tar_cells)
            if(length(tar_cells) > 1){
                dist = (1:ceiling((n_cells/2)))*1.8/n_cells + .25
                dist = c(dist, (1:(floor(n_cells/2)))*-1.8/n_cells + 1.15)
            } else {
                dist = .25
            }
            names(dist) = tar_cells
            tar_cells = tar_cells[which(dist > alphas[tar_cells])]
            if(length(tar_cells) != 0){
                alphas[tar_cells] = dist[tar_cells]
            }
            curve = curves[[j]]$s[names(times),]
            curve = curve[-c(1:400),]
            plot3d(curve[,1], curve[,2], curve[,3], type = 'l', add = TRUE, lwd = 2, 
                alpha = .5, col = 'black')
        }
        plot3d(p3[,1], p3[,2], p3[,3], box = FALSE, axes = FALSE, xlab = "", ylab = "", zlab = "", 
            col = cell_cols, alpha = alphas, add = TRUE, size = 1.5)
        bg3d(col = bg_color)
        par3d(windowRect = c(200, 300, 800, 800))
        view3d(userMatrix = rotationMatrix(2*pi*i/1000, 0, 1, 0))
        Sys.sleep(.1)
        rgl.snapshot(filename = paste0(out_dir,'/wflashing/', 
            sprintf("%04d", i), '.png'))
        rgl.close()
    
    }

    # Just spinning
    plot3d(p3[,1], p3[,2], p3[,3], box = FALSE, axes = FALSE, xlab = "", ylab = "", zlab = "", col = cell_cols, size = 1.5)
    par3d(windowRect = c(200, 300, 800, 800))
    bg3d(color = bg_color)

    for(i in 1:frame){
        view3d(userMatrix = rotationMatrix(2*pi*i/frame, 0, 1, 0))
        rgl.snapshot(filename = paste0(out_dir, '/wspinning/', 
            sprintf("%04d", i), '.png'))
    }
}