#' Aggregates the p-values from the two modalities (RNA and epigenomics) and then corrects them using the Benjamini-Hochberg method
#'
#' @param e.out1 a DataFrame (in the format of the emptyDrops outputs) with the following components:
#' Total:
#' Integer, the total count for each barcode.
#' LogProb:
#' Numeric, the log-probability of observing the barcode's count vector under the null model.
#' PValue:
#' Numeric, the Monte Carlo p-value against the null model.
#' Limited:
#' Logical, indicating whether a lower p-value could be obtained by increasing niters.
#' FDR:
#' Numeric, the p-values corrected using the Benjamini-Hochberg method
#'
#' @param e.out2 a DataFrame (in the format of the emptyDrops outputs) like e.out1
#'
#'
#' @return
#' a DataFrame (in the format of the emptyDrops outputs) with the following components:
#' Total_RNA:
#' Integer, the total count for each barcode.
#' LogProb_RNA:
#' Numeric, the log-probability of observing the barcode's count vector under the null model.
#' PValue_RNA:
#' Numeric, the Monte Carlo p-value against the null model.
#' Limited_RNA:
#' Logical, indicating whether a lower p-value could be obtained by increasing niters.
#' FDR_RNA:
#' Numeric, the p-values corrected using the Benjamini-Hochberg method
#' Total_chromatin:
#' Integer, the total count for each barcode.
#' LogProb_chromatin:
#' Numeric, the log-probability of observing the barcode's count vector under the null model.
#' PValue_chromatin:
#' Numeric, the Monte Carlo p-value against the null model.
#' Limited_chromatin:
#' Logical, indicating whether a lower p-value could be obtained by increasing niters.
#' FDR_chromatin:
#' Numeric, the p-values corrected using the Benjamini-Hochberg method
#' PValue_multi:
#' Numeric, the p-values of the droplets after aggregating PValue_chromatin and PValue_RNA for each droplet
#' FDR_multi:
#' Numeric, the result of the correction of PValue_multi using the Benjamini-Hochberg method
#'
#'
#' @examples
#'
#'
aggr_p_values <- function( e.out1, e.out2){
  # rename columns to before taking their disjoint union during merging
  colnames(e.out1) <- paste0(colnames(e.out1), "_RNA")
  colnames(e.out2) <- paste0(colnames(e.out2), "_chromatin")
  e_multi.out <- S4Vectors::merge(e.out1,e.out2, by=0, all=TRUE)  # by="col1"
  rownames(e_multi.out) <- e_multi.out$Row.names

  # aggregate p-values
  #aggr_pval <- apply( e_multi.out[,c('PValue_chromatin','PValue_RNA')], 1, function(x) fisher_v2(x) )
  aggr_pval <- apply( e_multi.out[,c('PValue_chromatin','PValue_RNA')], 1, function(x) balanced_aggr(x,r=1) )


  # correct aggregated p-values by the BH method to control the fdr
  fdr_of_aggr <- stats::p.adjust(aggr_pval, method="BH")

  e_multi.out["PValue_multi"] <- unname(aggr_pval)
  e_multi.out["FDR_multi"] <-  unname(fdr_of_aggr)

  return(e_multi.out)

}
