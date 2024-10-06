
#######################################################################
#######################################################################

#' Title Generate tidy/clean violin plot(s) from scRNA-seq analysis
#' @param seu_obj seurat object created with R package Seurat
#' @param features a vector containing genes
#' @param pt.size point size used for plotting
#' @param ncol arrange plots numbers in each row
#' @param cols colors used for plotting
#' @return combined_plot
#' @export
#' @examples
#' # Examples will be added in a future release.

# State dependencies and functions
#'@importFrom Seurat VlnPlot
#'@importFrom Seurat DefaultAssay
#'@importFrom ggplot2 labs
#'@importFrom ggplot2 theme
#'@importFrom patchwork wrap_plots

tidy.VlnPlot <- function(seu_obj, features, pt.size = 0.5, ncol = NULL,
                                  cols = NULL, split.by = NULL) {
  # Define data slot for analyzing
  DefaultAssay(seu_obj) <- "RNA"

  # Ensure features is a character vector
  features <- as.character(features)

  # Prepare a list to store plot(s)
  plots <- list()

  # Draw tidy violin plot(s)
  if (is.null(split.by)) {
  for (i in seq_along(features)) {
    plots[[i]] <- Seurat::VlnPlot(seu_obj, pt.size = pt.size, features = features[i], cols = cols)
    plots[[i]] <- plots[[i]] + labs(x = NULL) + labs(y = NULL) + NoLegend()
                                 }
         } else if (!is.null(split.by)) {
           for (i in seq_along(features)) {
                  plots[[i]] <- Seurat::VlnPlot(seu_obj, pt.size = pt.size, features = features[i], cols = cols,
                                               split.by = split.by, combine = TRUE)
                   plots[[i]] <- plots[[i]] + labs(x = NULL) + labs(y = NULL)
                                      }
                                                   }
  # Arrange the plot(s) layout
   if (is.null(ncol)) {
      if(length(features) <= 3){
         ncol = length(features)} else if (length(features) == 4 ) {
           ncol = 2} else if (length(features) == 5 | length(features) == 6) {
              ncol = 3} else if (length(features) == 7 | length(features) == 8) {
                 ncol = 4} else if(length(features) == 9) {
                    ncol = 3} else if (length(features) == 10) {
                      ncol = 5} else if (length(features) >= 11 && length(features) <= 14) {
                        ncol = 4} else if (length(features) == 15) {
                          ncol = 5} else if(length(features) == 16) {
                            ncol = 4} else if(length(features) >= 17) {
                              ncol = 6}
                                }

  # Combine plot(s) into one single sheet
  combined_plot <- wrap_plots(plots, ncol = ncol)

  # Remove redundant object to save memory
  rm(plots)

  # Return ploting results
  print(combined_plot)
     }
#######################################################################
#######################################################################





#######################################################################
#######################################################################
#' Title Generate tidy/clean feature plot(s) from scRNA-seq analysis
#' @param seu_obj seurat object created with R package Seurat
#' @param features a vector containing genes
#' @param pt.size point size used for plotting
#' @param ncol arrange plots numbers in each row
#' @param cols colors used for plotting
#' @return combined_plot
#' @export
#' @examples
#' # Examples will be added in a future release.

# State dependencies and functions
#'@importFrom Seurat FeaturePlot
#'@importFrom Seurat DefaultAssay
#'@importFrom ggplot2 labs
#'@importFrom ggplot2 theme
#'@importFrom patchwork wrap_plots

tidy.FeaturePlot <- function(seu_obj, features, pt.size = 0.5, ncol = NULL,
                             cols = c("gray", "red"), Legend = TRUE) {
  # Define data slot for analyzing
  DefaultAssay(seu_obj) <- "RNA"

  # Ensure features is a character vector
  features <- as.character(features)

  # Prepare a list to store plot(s)
  plots <- list()

  # Draw tidy feature plot(s)
  for (i in seq_along(features)) {
    plots[[i]] <- Seurat::FeaturePlot(seu_obj, pt.size = pt.size, features = features[i], cols = cols) +
      NoAxes()  +
      theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
                                 }
  # Handle legends in the plot(s)
  if (isTRUE(Legend)) {
    for (i in seq_along(features)) {
             plots[[i]] <- plots[[i]]
                                  }
                     } else if (isFALSE(Legend)) {
                          for (i in seq_along(features)) {
                            plots[[i]] <- plots[[i]] + NoLegend()
                          }
                     }


  # Arrange the plot(s) layout
  if (is.null(ncol)) {
    if(length(features) <= 3){
      ncol = length(features)} else if (length(features) == 4 ) {
        ncol = 2} else if (length(features) == 5 | length(features) == 6) {
          ncol = 3} else if (length(features) == 7 | length(features) == 8) {
            ncol = 4} else if(length(features) == 9) {
              ncol = 3} else if (length(features) == 10) {
                ncol = 5} else if (length(features) >= 11 && length(features) <= 14) {
                  ncol = 4} else if (length(features) == 15) {
                    ncol = 5} else if(length(features) == 16) {
                      ncol = 4} else if(length(features) >= 17) {
                        ncol = 6}
                      }

  # Combine plot(s) into one single sheet
  combined_plot <- wrap_plots(plots, ncol = ncol)

  # Remove redundant object to save memory
  rm(plots)

  # Return ploting results
  print(combined_plot)
      }
#######################################################################
#######################################################################






