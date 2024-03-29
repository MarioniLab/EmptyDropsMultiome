% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emptydrops_multiome.R
\name{emptydrops_multiome}
\alias{emptydrops_multiome}
\title{Calls cells using emptydrops in two genomic modalities.}
\usage{
emptydrops_multiome(
  count_matrix_rna,
  lower_rna = NULL,
  barhop_rna = NULL,
  count_matrix_atac,
  lower_atac = NULL,
  barhop_atac = NULL,
  niter_rna = 10000,
  niter_atac = 25000,
  verbose = TRUE
)
}
\arguments{
\item{count_matrix_rna}{dgCMatrix with the counts for the rna modality.}

\item{lower_rna}{single real number at or below which droplets are used to model the RNA soup.}

\item{barhop_rna}{single real number below which droplets are assumed to be full of barhops and they are excluded from modeling the soup.}

\item{count_matrix_atac}{dgCMatrix with the counts for the ATAC modality.}

\item{lower_atac}{single real number at or below which droplets are used to model the ATAC soup.}

\item{barhop_atac}{single real number below which droplets are assumed to be full of barhops and they are excluded from modeling the soup.}

\item{niter_rna}{single real number for the number of iteration to be performed to statistically compare the RNA modality with the ambient profile.}

\item{niter_atac}{single real number for the number of iteration to be performed to statistically compare the ATAC modality with the ambient profile.}

\item{verbose}{if TRUE various intermediate steps and plots are printed.}
}
\value{
a DataFrame object (from the S4Vectors package) with the following components:
Total_RNA:
Integer, the total count for each barcode.
LogProb_RNA:
Numeric, the log-probability of observing the barcode's count vector under the null model.
PValue_RNA:
Numeric, the Monte Carlo p-value against the null model.
Limited_RNA:
Logical, indicating whether a lower p-value could be obtained by increasing niters.
FDR_RNA:
Numeric, the p-values corrected using the Benjamini-Hochberg method
Total_chromatin:
Integer, the total count for each barcode.
LogProb_chromatin:
Numeric, the log-probability of observing the barcode's count vector under the null model.
PValue_chromatin:
Numeric, the Monte Carlo p-value against the null model.
Limited_chromatin:
Logical, indicating whether a lower p-value could be obtained by increasing niters.
FDR_chromatin:
Numeric, the p-values corrected using the Benjamini-Hochberg method
PValue_multi:
Numeric, the p-values of the droplets after aggregating PValue_chromatin and PValue_RNA for each droplet
FDR_multi:
Numeric, the result of the correction of PValue_multi using the Benjamini-Hochberg method
k_means:
Logical, indicating whether droplet lies above the k-means line
above_ambiguous:
Logical, indicating whether droplet lies above the ambiguous area, where droplets are assumed to unambiguously contain cells.
}
\description{
Calls cells using emptydrops in two genomic modalities.
}
