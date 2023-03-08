#' Calls cells using emptydrops on a single modality. It models the soup using droplets with counts between barhop_end and lower.
#'
#' @param count_matrix dgCMatrix with the counts.
#' @param lower single real number at or below which droplets are used to model the RNA soup.
#' @param barhop_end single real number below which droplets are assumed to be full of barhops and they are excluded from modeling the soup.
#' @param niters the number of iterations to use for the simulations that compare the profile of each droplet to the ambient profile.
#'
#' @return a DataFrame (in the format of the emptyDrops outputs) with the following components:
#' Total*:
#' Integer, the total count for each barcode.
#' LogProb:
#' Numeric, the log-probability of observing the barcode's count vector under the null model.
#' PValue:
#' Numeric, the Monte Carlo p-value against the null model.
#' Limited:
#' Logical, indicating whether a lower p-value could be obtained by increasing niters.
#' FDR:
#' Numeric, the p-values corrected using the Benjamini-Hochberg method
#' @export
#'
#' @examples
emptydrops_with_barhop <- function(count_matrix, lower, barhop_end, niters=10000){
  # calls cells using a single modality
  # inputs: count_matrix is a features-by-cell matrix,
  #         lower is the upper bound on the count for a droplet to be used to model the soup
  #         barhop_end is the lower bound on the count for a droplet to be used to model the soup
  # outputs: dataframe like e.out


  # # inputs by hand (to be automated)
  # barhop_end = 20
  # lower = 120

  # subset out the barhops from count matrix to feed it into emptydrops
  #barhops <- colnames(count_matrix)[colSums(count_matrix)<barhop_end]
  #counts_wo_barhop <- count_matrix[,! colnames(count_matrix) %in% barhops]
  barhops <- Matrix::colSums(count_matrix)<barhop_end
  counts_wo_barhop <- count_matrix[, !barhops]
  if (!length(colnames(count_matrix))==0){
    barhop_names <- colnames(count_matrix)[barhops]
  }else{
    barhop_names <- NULL
  }

  # run emptyrdops on counts_wo_barhop
  e.out <- DropletUtils::emptyDrops(counts_wo_barhop, lower=lower, niters=niters, test.ambient=TRUE  )

  # create df with NA for barhops
  totals_barhop <- Matrix::colSums(count_matrix)[barhops]
  totals_barhop <- unname(totals_barhop)
  all.lr_barhops <- rep(as.numeric(NA), sum(barhops))
  all.p_barhops <- rep(as.numeric(NA), sum(barhops))
  all.lim_barhops <- rep(as.numeric(NA), sum(barhops))
  #astats$metadata <- rep(as.numeric(NA), length(barhops))
  fdr_barhop <- rep(as.numeric(NA), sum(barhops))

  # output <- DataFrame(Total=totals, LogProb=all.lr, PValue=all.p, Limited=all.lim, row.names=colnames(m))
  # metadata(output) <- c(astats$metadata, list(niters=niters, ambient=ambient.prop, alpha=alpha))
  #barhop_output <- DataFrame(Total=totals_barhop, LogProb=all.lr_barhops, PValue=all.p_barhops, Limited=all.lim_barhops, row.names=barhops, FDR=fdr_barhop)
  barhop_output <- S4Vectors::DataFrame(Total=totals_barhop, LogProb=all.lr_barhops, PValue=all.p_barhops, Limited=all.lim_barhops, row.names=barhop_names, FDR=fdr_barhop)
  #metadata(barhop_output) <- c(astats$metadata, list(niters=niters, ambient=ambient.prop, alpha=alpha))

  # merge the two dfs
  eD.out <- rbind(barhop_output, e.out)

  eD.out["barhop_amb_tent"] <- as.numeric(eD.out$Total<barhop_end) + as.numeric(eD.out$Total<=lower) + as.numeric(eD.out$Total>lower) *3

  return(eD.out)
}

