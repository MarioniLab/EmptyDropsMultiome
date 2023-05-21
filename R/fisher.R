#' Use chi-squared function to aggregate p-values derived from different hypothesis tests
#'
#' @param pvalues Vector of real numbers.
#'
#' @return Real number of the aggregated p-value.
#'
fisher_v2 <- function(pvalues)
{
  # use chi-squared function to aggregate p-values derived from different hypothesis tests
  # input: pvalues vector of real numbers
  # output: real number of the aggregated p-value
  if (any(is.na(pvalues)) )
  {
    return(NA)
  }
  if(length(pvalues)==0)
  {
    return(NA)
  }
  if(any(pvalues <0) || any(pvalues >1))
  {
    stop('p-values must be in [0, 1]')
  }
  if(length(pvalues)==1)
  {
    return(pvalues)
  }

  chisq = -2 * sum(log(pvalues))
  dof <- 2* length(pvalues)
  stats::pchisq(chisq, dof, lower.tail = FALSE)
}



#' Mean aggregation of real values.
#'
#' @param pvalues Vector of real numbers.
#' @param r Real number for the kind of mean which will be used to aggregate the p-values.
#' For instance, r=1 is the arithmetic mean, r=-1 is the harmonic mean etc.
#'
#' @return Real number of the aggregated p-value.
#'
mean_aggr <- function(pvalues, r){
  if(r==0)
  {
    stop('r must be non zero')
  }
  if (is.na(r) )
  {
    stop('r must have a value')
  }
  if (any(is.na(pvalues)) )
  {
    return(NA)
  }

  aggr <- (  sum(pvalues^r)/length(pvalues)  )^(1/r)
  return(aggr)

}

#' Aggregate p-values using the average of the Fisher aggregation value and the mean aggregation method.
#'
#' @param pvalues Vector of real numbers in the range from 0 to 1.
#' @param r r Real number for the kind of mean which will be used to aggregate the p-values.
#' For instance, r=1 is the arithmetic mean, r=-1 is the harmonic mean etc.
#'
#' @return Real number of the aggregated p-value.
#'
balanced_aggr <- function(pvalues, r){
  if(r==0)
  {
    stop('r must be non zero')
  }
  if (is.na(r) )
  {
    stop('r must have a value')
  }
  if (any(is.na(pvalues)) )
  {
    return(NA)
  }

  #balanced <- (  mean_aggr(pvalues, r) + fisher_v2(pvalues)  ) /2
  balanced <- mean_aggr(pvalues, r)
  return(balanced)

}
