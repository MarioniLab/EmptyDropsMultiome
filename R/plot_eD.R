#' Creates the cell calling plot.
#'
#' @param x_axis Vector with values for the x-coordinate.
#' @param y_axis Vector with values for the y-coordinate.
#' @param colors Vector with values for the color to be used in the scatter plot.
#' @param title character to be used as title for the plot.
#'
#' @return invisible(NULL)
#' @export
#'
#' @examples
cell_calling_plot <- function(x_axis, y_axis, colors, title){


  plot(x_axis, y_axis,
       #     xlim=c(0,1510) , ylim=c(0,8),
       pch=".",
       # cex=2,
       col=colors,
       xlab="log10(atac_count)", ylab="log10(rna_count)",
       main=paste0(title)
  )+ggplot2::theme(axis.title.x=ggplot2::element_text(size=14,face="bold"),
          axis.title.y=ggplot2::element_text(size=14,face="bold"))


  return(invisible(NULL))

}
