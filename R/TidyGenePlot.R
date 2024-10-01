
########################################################
########################################################
# generate tidy violin plot(s) from scRNA-seq analysis
# seu_obj:seurat object created by the R package Seurat from scRNA-seq analysis
# fatures: a vector of genes(not a list)
# pt.size: point size for plots (must be between 0 and 1)
# ncol: number of columns to arrange plots in each row


#' Title
#'
#' @param seu_obj seurat object created with R package Seurat
#' @param features a vector containing genes
#' @param pt.size point size used for plotting
#' @param ncol arrange plots numbers in each row
#' @param cols colors used for plotting
#'
#' @return combined_plot
#' @export
#' @examples
#'@importFrom Seurat VlnPlot
#'@importFrom Seurat FeaturePlot
#'@importFrom patchwork wrap_plots
tidy.VlnPlot <- function(seu_obj, features, pt.size = 0.5, ncol = 1, cols = NULL) {

  DefaultAssay(seu_obj) <- "RNA"

  # ensure features is a character vector
  features <- as.character(features)

  # prepare a list to store plot(s)
  plots <- list()

  # draw tidy violin plot(s)

  for (i in seq_along(features)) {
    plots[[i]] <- Seurat::VlnPlot(seu_obj, pt.size = pt.size, features = features[i], cols = cols)
    plots[[i]] <- plots[[i]] + labs(x = NULL) + labs(y = NULL) + NoLegend()
  }
  combined_plot <- wrap_plots(plots, ncol = ncol)
  rm(plots)
  print(combined_plot)

}




#' Title
#'
#' @param seu_obj seurat object created with R package Seurat
#' @param features a vector containing genes
#' @param pt.size point size used for plotting
#' @param ncol arrange plots numbers in each row
#' @param cols colors used for plotting
#'
#' @return combined_plot
#' @export
#'
#' @examples
#'
tidy.FeaturePlot <- function(seu_obj, features, pt.size = 0.5, ncol = 1,
                             cols = c("gray", "blue")) {

  DefaultAssay(seu_obj) <- "RNA"

  # Ensure features is a character vector
  features <- as.character(features)

  # Prepare a list to store plot(s)
  plots <- list()

  # Draw tidy feature plot(s)
  for (i in seq_along(features)) {
    plots[[i]] <- Seurat::FeaturePlot(seu_obj, pt.size = pt.size, features = features[i], cols = cols) + NoAxes() +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  }
  combined_plot <- wrap_plots(plots, ncol = ncol)
  rm(plots)
  print(combined_plot)
}







