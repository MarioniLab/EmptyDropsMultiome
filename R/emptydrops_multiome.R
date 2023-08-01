#' Calls cells using emptydrops in two genomic modalities.
#'
#' @param count_matrix_rna dgCMatrix with the counts for the rna modality.
#' @param lower_rna single real number at or below which droplets are used to model the RNA soup.
#' @param barhop_rna single real number below which droplets are assumed to be full of barhops and they are excluded from modeling the soup.
#' @param count_matrix_atac dgCMatrix with the counts for the ATAC modality.
#' @param lower_atac single real number at or below which droplets are used to model the ATAC soup.
#' @param barhop_atac single real number below which droplets are assumed to be full of barhops and they are excluded from modeling the soup.
#' @param niter_rna single real number for the number of iteration to be performed to statistically compare the RNA modality with the ambient profile.
#' @param niter_atac single real number for the number of iteration to be performed to statistically compare the ATAC modality with the ambient profile.
#' @param verbose if TRUE various intermediate steps and plots are printed.
#'
#'
#'
#' @return a DataFrame object (from the S4Vectors package) with the following components:
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
#' k_means:
#' Logical, indicating whether droplet lies above the k-means line
#' above_ambiguous:
#' Logical, indicating whether droplet lies above the ambiguous area, where droplets are assumed to unambiguously contain cells.
#'
#'
#' @export
#'
emptydrops_multiome <- function(count_matrix_rna, lower_rna=NULL, barhop_rna=NULL, count_matrix_atac, lower_atac=NULL, barhop_atac=NULL, niter_rna=10000, niter_atac=25000, verbose=TRUE, seed=42 ){
  
  if (!is.null(lower_rna) & !is.null(barhop_rna) & !is.null(lower_atac) & !is.null(barhop_atac)  ){
    if (lower_rna<0 | barhop_rna<0 | lower_atac<0 | barhop_atac<0){
      stop("barhop and lower parameters should be non negative")
    }
  }

  if ( is.null(lower_rna)  | is.null(barhop_rna) ){
    exp_counts_per_cell = unname(Matrix::colSums(count_matrix_rna))
    mu1 = fit_3_normals(exp_counts_per_cell, verbose)
    lower_rna = mu1[1]
    barhop_rna = mu1[2]
    print(paste0("the rna lower is ",lower_rna ) )
    print(paste0("the rna barhop is ",barhop_rna ) )
  }else{
    exp_counts_per_cell_rna = unname(Matrix::colSums(count_matrix_rna))
    observations = data.frame( "nCount_RNA" = exp_counts_per_cell_rna )

    print(ggplot2::ggplot(observations, ggplot2::aes(x = nCount_RNA)) +
      ggplot2::geom_histogram(binwidth =1) +
      ggplot2::xlim(0, 700)+
      ggplot2::scale_y_continuous( limit = c(0, 100000), oob = function(x, limits) x)+
      ggplot2::geom_vline(xintercept = lower_rna,
                          # linetype="dotted",
                          color = "blue",
                          size=0.5)+
      ggplot2::geom_vline(xintercept = barhop_rna,
                          # linetype="dotted",
                          color = "red",
                          size=0.5) )


  }

  if (  is.null(lower_atac) | is.null(barhop_atac)  ){
    exp_counts_per_cell = unname(Matrix::colSums(count_matrix_atac))
    
    mu2 = fit_3_normals_atac( exp_counts_per_cell, verbose )
    lower_atac = mu2[1]
    barhop_atac = mu2[2]
    print(paste0("the atac lower is ", lower_atac ) )
    print(paste0("the atac barhop is ", barhop_atac ) )
  }else{
    exp_counts_per_cell = unname(Matrix::colSums(count_matrix_atac))
    observations = data.frame( "nCount_ATAC" = exp_counts_per_cell )

    ggplot2::ggplot(observations, ggplot2::aes(x = nCount_ATAC)) +
      ggplot2::geom_histogram(binwidth =1) +
      ggplot2::scale_x_continuous( limit = c(0, 500), oob = function(x, limits) x)+
      ggplot2::scale_y_continuous( limit = c(0, 100000), oob = function(x, limits) x)+
      ggplot2::geom_vline(xintercept = lower_atac,
                          # linetype="dotted",
                          color = "blue",
                          size=0.5)+
      ggplot2::geom_vline(xintercept = barhop_atac,
                          # linetype="dotted",
                          color = "red",
                          size=0.5)


  }

  # call cells based on RNA
  eD.out_rna <- emptydrops_with_barhop(count_matrix_rna, lower_rna, barhop_rna, niters=niter_rna, seed=seed)

  # call cells based on epigenetics
  eD.out_atac <- emptydrops_with_barhop(count_matrix_atac, lower_atac, barhop_atac, niters=niter_atac, seed=seed)

  # call cells based on both RNA and epigenetics
  e_multi.out <- aggr_p_values( eD.out_rna, eD.out_atac)

  # accept cells based on k_means
  e_multi.out <- accept_k_means(e_multi.out, lower_atac, lower_rna)

  # add metadata to e_multi.out
  # e_multi.out$lower_rna = lower_rna
  # e_multi.out$barhop_rna=barhop_rna
  # e_multi.out$lower_atac =lower_atac
  # e_multi.out$barhop_atac=barhop_atac

  e_multi.out@metadata["lower_rna"] = lower_rna
  e_multi.out@metadata["barhop_rna"]=barhop_rna
  e_multi.out@metadata["lower_atac"] =lower_atac
  e_multi.out@metadata["barhop_atac"]=barhop_atac

  return(e_multi.out)

}


#' #' Plots histogram of counts
#' #'
#' #' @param ct_matrix: a matrix of counts
#' #' @param assay_name: the type of counts e.g. "RNA" or "ATAC"
#' #' @param lower
#' #' @param barhop_end
#' #'
#' #' @return plots
#' #'
#' ct_hist <- function(ct_matrix, assay_name, lower, barhop_end){
#'
#'   exp_counts_per_cell = unname(Matrix::colSums(ct_matrix))
#'   observations = data.frame( "nCount" = exp_counts_per_cell )
#'   colnames(observations) <- paste0("nCount_", assay_name)
#'
#'   print(ggplot2::ggplot(observations, ggplot2::aes(x = nCount_RNA)) +
#'           ggplot2::geom_histogram(binwidth =1) +
#'           ggplot2::xlim(0, 500)+
#'           ggplot2::scale_y_continuous( limit = c(0, 80000), oob = function(x, limits) x)+
#'           ggplot2::geom_vline(xintercept = lower_rna,
#'                               # linetype="dotted",
#'                               color = "blue",
#'                               size=0.5)+
#'           ggplot2::geom_vline(xintercept = barhop_rna,
#'                               # linetype="dotted",
#'                               color = "red",
#'                               size=0.5) )
#'
#'   ggplot2::ggplot(observations, ggplot2::aes(x = nCount_RNA)) +
#'     ggplot2::geom_histogram(binwidth =1) +
#'     ggplot2::scale_x_continuous( limit = c(0, 10000), oob = function(x, limits) x)+
#'     ggplot2::scale_y_continuous( limit = c(0, 40), oob = function(x, limits) x)+
#'     ggplot2::geom_vline(xintercept = lower_rna,
#'                         # linetype="dotted",
#'                         color = "blue",
#'                         size=0.5)+
#'     ggplot2::geom_vline(xintercept = barhop_rna,
#'                         # linetype="dotted",
#'                         color = "red",
#'                         size=0.5)
#'
#'
#'
#' }

