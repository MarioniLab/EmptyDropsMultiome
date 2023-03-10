#' After using emptydrops, set FDR_multi<-1 for cells above the ambiguous area, and set FDR_multi<-0 in cells below the lower line.
#'
#' @param e_multi.out Dataframe with columns for Total_RNA, Total_chromatin, rownames, barhop_amb_tent_chromatin, barhop_amb_tent_rna, FDR_multi
#' @param lower_atac A single number to act as the value of the lower parameter in ATAC
#' @param lower_rna A single number to act as the value of the lower parameter in RNA
#'
#' @return Dataframe with columns for Total_RNA, Total_chromatin, rownames, barhop_amb_tent_chromatin,
#'         barhop_amb_tent_rna, FDR_multi, FDR, above_ambiguous, k_means
#' @export
#'
accept_k_means <- function( e_multi.out, lower_atac, lower_rna){
  # it changes to 0 the FDR_multi value of droplets that lie above the k-means threshold
  # input: e_multi.out with a column "FDR_multi"
  # output: e_multi.out with adjusted column "FDR_multi"

  # create count dataframe with columns 'atac_count' (int), 'rna_count' (int), 'excluded' (bool), 'is_cell'

  count_df <- data.frame("rownames" = rownames(e_multi.out),
                         "atac_count" = e_multi.out$Total_chromatin,
                         "rna_count" = e_multi.out$Total_RNA,
                         "excluded"  = rep(F, length(e_multi.out$Total_RNA)),
                         "is_cell"  =  rep(0, length(e_multi.out$Total_RNA))
  )
  count_df$atac_count[is.na(count_df$atac_count)] <- 0
  count_df$rna_count[is.na(count_df$rna_count)] <- 0

  list_output <- call_cells(count_df)
  count_df_new <- list_output[[1]]
  equation_k_means <- list_output[[2]]

  e_multi.out_new <- e_multi.out

  # replace NA values with 0
  e_multi.out_new$Total_RNA[is.na(e_multi.out_new$Total_RNA)] = 0
  e_multi.out_new$Total_chromatin[is.na(e_multi.out_new$Total_chromatin)] = 0

  # tag according to k_means (ie cellranger-arc)
  cell_barcodes <- count_df_new$rownames[which(count_df_new$is_cell==1) ]
  e_multi.out_new$k_means <- 0
  e_multi.out_new[cell_barcodes, ]$k_means <- 1

  # calculate the ambiguous_above area and tag cells
  equation_parallel <- calc_ambiguous_above(e_multi.out_new[,c("Total_RNA", "Total_chromatin")], equation_k_means)
  # cells_above_ambiguous <- e_multi.out_new$Total_chromatin - equation_k_means[2] * e_multi.out_new$Total_RNA > equation_k_means[1] & e_multi.out_new$Total_chromatin - equation_k_means[2] * e_multi.out_new$Total_RNA > equation_parallel[1]
  print(equation_parallel)
  cells_above_ambiguous <-  log10(e_multi.out_new$Total_RNA+0.1) - equation_parallel[2] * log10(e_multi.out_new$Total_chromatin + 0.1) > equation_parallel[1]
  e_multi.out_new$above_ambiguous <- 0
  e_multi.out_new[cells_above_ambiguous, ]$above_ambiguous <- 1



  # cell calling plots
  graphics::par(mfrow = c(2, 2)) # Create a 2 x 2 plotting matrix
  graphics::par(mar=c(5,5,4,1)+.1)

  cell_calling_plot(log10(e_multi.out$Total_chromatin + 0.1), log10(e_multi.out$Total_RNA + 0.1), factor(e_multi.out$FDR_multi<0.001 ), "eD_multi tentative")
  graphics::abline(a=equation_k_means[1],b=equation_k_means[2], col="blue")
  graphics::abline(a=equation_parallel[1],b=equation_parallel[2], col="blue")

  cell_calling_plot(log10(e_multi.out$Total_chromatin + 0.1), log10(e_multi.out$Total_RNA + 0.1), factor(rownames(e_multi.out) %in% cell_barcodes ), "cR_multi calling")
  graphics::abline(a=equation_k_means[1],b=equation_k_means[2], col="blue")

  # reject by default droplets below both lowers or at least one of the barhops or the na
  e_multi.out_new$FDR_multi[e_multi.out_new$barhop_amb_tent_RNA==2] <- 1
  e_multi.out_new$FDR_multi[e_multi.out_new$barhop_amb_tent_chromatin==2] <- 1


  # calculate equation_lower and reject cells below it
  equation_lower = calc_lower_line(log10(lower_atac+0.1), log10(lower_rna+0.1), equation_parallel)
  e_multi.out_new$FDR_multi[log10(e_multi.out_new$Total_RNA+0.1) < equation_lower[2] * log10(e_multi.out_new$Total_chromatin + 0.1) + equation_lower[1] ] <- 1

  # turn FDR_multi==NA to FDR_multi==1
  e_multi.out_new$FDR_multi[is.na(e_multi.out_new$FDR_multi)] <- 1

  # plot effect of rejecting below lower and FDR_multi NA->1
  cell_calling_plot(log10(e_multi.out_new$Total_chromatin + 0.1), log10(e_multi.out_new$Total_RNA + 0.1), factor(e_multi.out_new$FDR_multi<0.001 ), "eD_multi after rejecting below lower")
  graphics::abline(a=equation_k_means[1],b=equation_k_means[2], col="blue")
  graphics::abline(a=equation_parallel[1],b=equation_parallel[2], col="blue")
  graphics::abline(a=equation_lower[1],b=equation_lower[2], col="blue")


  # accept by default all the droplets above the ambiguous area
  e_multi.out_new$FDR <- e_multi.out_new$FDR_multi
  e_multi.out_new[cells_above_ambiguous, ]$FDR_multi <- 0

  cell_calling_plot(log10(e_multi.out_new$Total_chromatin + 0.1), log10(e_multi.out_new$Total_RNA + 0.1), factor(e_multi.out_new$FDR_multi<0.001 ), "eD_multi, after accepting above knee (final)")
  graphics::abline(a=equation_k_means[1],b=equation_k_means[2], col="blue")
  graphics::abline(a=equation_parallel[1],b=equation_parallel[2], col="blue")
  graphics::abline(a=equation_lower[1],b=equation_lower[2], col="blue")

  e_multi.out_new@metadata["k_means_slope"] = equation_k_means[2]
  e_multi.out_new@metadata["k_means_intercept"] = equation_k_means[1]

  return(e_multi.out_new)
}
