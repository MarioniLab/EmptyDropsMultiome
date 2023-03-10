#' Fits a mixture of three normal distributions to histogram of RNA counts, to infer the parameters lower and barhop_end
#' which are used in emptydrops.
#'
#' @param exp_counts_per_cell a vector with RNA library size per barcode.
#' @param verbose a Boolean variable indicating whether to print computation messages.
#'
#' @return A 2d vector with the values lower and barhop_end for RNA.
#'
#' @export
#'
#' @examples
fit_3_normals <- function(exp_counts_per_cell, verbose){

  #exp_counts_per_cell = unname(Matrix::colSums(srat_029_all@assays[["RNA"]]@counts))
  #exp_counts_per_cell = unname(Matrix::colSums(count_matrix))

  my_mix <- mixtools::normalmixEM(exp_counts_per_cell, mu=c(5,500,1200), sigma=c(5,50,800) , epsilon=1e-09)

  mu_order = order(my_mix$mu)
  mu_barhop = my_mix$mu[mu_order[1]]
  mu_ambient = my_mix$mu[mu_order[2]]
  mu_true_cells = my_mix$mu[mu_order[3]]
  lower = mu_ambient+1.5*my_mix$sigma[mu_order[2]]
  barhop = mu_barhop+4*my_mix$sigma[mu_order[1]]
  
  n = length(exp_counts_per_cell)
  equil_pt = equil_of_normals(mu_barhop, my_mix$sigma[mu_order[1]], my_mix$lambda[mu_order[1]],
                              mu_ambient, my_mix$sigma[mu_order[2]], my_mix$lambda[mu_order[2]])
  
  if (verbose){
    
    print(paste0("the mu's are : ", my_mix$mu[1],", ", my_mix$mu[2],", ", my_mix$mu[3] ))
    my_mix$mu
    print(paste0("the sigma's are : ", my_mix$sigma[1],", ", my_mix$sigma[2],", ", my_mix$sigma[3] ))
    my_mix$sigma
    print(paste0("the lambdas are : ", my_mix$lambda[1],", ", my_mix$lambda[2],", ", my_mix$lambda[3] ))
    my_mix$lambda
    print(paste0("equil point is ", equil_pt) )
    print(paste0("the lower is : ", lower ))
    print(paste0("the barhop_end is : ", barhop ))
    
  }
  
  observations = data.frame( "nCount_RNA" = exp_counts_per_cell )

  q1 = ggplot2::ggplot(observations, ggplot2::aes(x = nCount_RNA)) +
    ggplot2::geom_histogram(binwidth =1) +
    mapply(
      function(mean, sd, lambda, n, binwidth) {
        ggplot2::stat_function(
          fun = function(x) {
            (stats::dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
          }
        )
      },
      mean = my_mix[["mu"]], #mean
      sd = my_mix[["sigma"]], #standard deviation
      lambda = my_mix[["lambda"]], #amplitude
      n = length(observations$nCount_RNA), #sample size
      binwidth = 1  #binwidth used for histogram
    ) +
    ggplot2::scale_x_continuous( limit = c(0, 800), oob = function(x, limits) x)+
    ggplot2::scale_y_continuous( limit = c(0, 15000), oob = function(x, limits) x)+
    ggplot2::geom_vline(xintercept = lower,
               # linetype="dotted",
               color = "blue",
               size=0.5)+
    ggplot2::geom_vline(xintercept = barhop,
               # linetype="dotted",
               color = "red",
               size=0.5)

  q2 = ggplot2::ggplot(observations, ggplot2::aes(x = nCount_RNA)) +
    ggplot2::geom_histogram(binwidth =1) +
    mapply(
      function(mean, sd, lambda, n, binwidth) {
        ggplot2::stat_function(
          fun = function(x) {
            (stats::dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
          }
        )
      },
      mean = my_mix[["mu"]], #mean
      sd = my_mix[["sigma"]], #standard deviation
      lambda = my_mix[["lambda"]], #amplitude
      n = length(observations$nCount_RNA), #sample size
      binwidth = 1  #binwidth used for histogram
    ) +
    ggplot2::scale_x_continuous( limit = c(0, lower+30), oob = function(x, limits) x)+
    ggplot2::scale_y_continuous( limit = c(0, 10000), oob = function(x, limits) x)+
    ggplot2::geom_vline(xintercept = lower,
               # linetype="dotted",
               color = "blue",
               size=0.5)+
    ggplot2::geom_vline(xintercept = barhop,
               # linetype="dotted",
               color = "red",
               size=0.5)


  q3 = ggplot2::ggplot(observations, ggplot2::aes(x = nCount_RNA)) +
    ggplot2::geom_histogram(binwidth =1) +
    mapply(
      function(mean, sd, lambda, n, binwidth) {
        ggplot2::stat_function(
          fun = function(x) {
            (stats::dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
          }
        )
      },
      mean = my_mix[["mu"]], #mean
      sd = my_mix[["sigma"]], #standard deviation
      lambda = my_mix[["lambda"]], #amplitude
      n = length(observations$nCount_RNA), #sample size
      binwidth = 1  #binwidth used for histogram
    ) +
    ggplot2::scale_x_continuous( limit = c(true_cells-my_mix$sigma[mu_order[3]], true_cells+my_mix$sigma[mu_order[3]]), oob = function(x, limits) x)+
    ggplot2::ylim(0,75)

  egg::ggarrange(q1, q2, q3 ,
                     labels = c("A", "B", "C"),
                     ncol = 3, nrow = 1)

  return(c(lower, barhop))
}







#' Fits a mixture of three normal distributions to histogram of ATAC counts, to infer the parameters lower and barhop_end
#' which are used in emptydrops.
#'
#' @param exp_counts_per_cell a vector with the ATAC library size per barcode.
#'
#' @return A 2d vector with the values lower and barhop_end for ATAC.
#'
#' @export
#'
#' @examples
fit_3_normals_atac <- function(exp_counts_per_cell){
  
  #exp_counts_per_cell = unname(Matrix::colSums(count_matrix))
  
  
  #my_mix <- normalmixEM(exp_counts_per_cell, mu=c(5,500,1200), sigma=c(5,50,800) , epsilon=1e-09)
  t <- try(my_mix <- mixtools::normalmixEM(exp_counts_per_cell, lambda=c(6*10^5/(7*10^5), 8*10^4/(7*10^5), 10^4/(7*10^5) ), 
                                           mu=c(5,100,20000), sigma=c(5,50,800) , epsilon=1e-09))
  print(inherits(t,"try-error"))
  if( inherits(t,"try-error")  ){
    my_mix <- data.frame("mu"=c(1, 40, 800  ), "sigma"=c(1, 60, 600  ), "lambda"= c(1,0.10,0.05))
  }

  print(paste0("the mus are : ", my_mix$mu[1],", ", my_mix$mu[2],", ", my_mix$mu[3] ))
  my_mix$mu
  print(paste0("the sigmas are : ", my_mix$sigma[1],", ", my_mix$sigma[2],", ", my_mix$sigma[3] ))
  my_mix$sigma
  print(paste0("the lambdas are : ", my_mix$lambda[1],", ", my_mix$lambda[2],", ", my_mix$lambda[3] ))
  my_mix$lambda
  mu_order = order(my_mix$mu)
  mu_order

  lower = my_mix$mu[mu_order[2]]+2*my_mix$sigma[mu_order[2]]
  if (lower>180){
    lower = 180
  }
  lower
  barhop = my_mix$mu[mu_order[1]]+4*my_mix$sigma[mu_order[1]]
  barhop
  true_cells = my_mix$mu[mu_order[3]]
  observations = data.frame( "nCount_ATAC" = exp_counts_per_cell )

  q1 =ggplot2::ggplot(observations, ggplot2::aes(x = nCount_ATAC)) +
    ggplot2::geom_histogram(binwidth =1) +
    mapply(
      function(mean, sd, lambda, n, binwidth) {
        ggplot2::stat_function(
          fun = function(x) {
            (stats::dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
          }
        )
      },
      mean = my_mix[["mu"]], #mean
      sd = my_mix[["sigma"]], #standard deviation
      lambda = my_mix[["lambda"]], #amplitude
      n = length(observations$nCount_ATAC), #sample size
      binwidth = 1  #binwidth used for histogram
    ) +
    ggplot2::scale_x_continuous( limit = c(0, 500), oob = function(x, limits) x)+
    ggplot2::scale_y_continuous( limit = c(0, 15000), oob = function(x, limits) x)+
    ggplot2::geom_vline(xintercept = lower,
               # linetype="dotted",
               color = "blue",
               size=0.5)+
    ggplot2::geom_vline(xintercept = barhop,
               # linetype="dotted",
               color = "red",
               size=0.5)
  
  n = length(observations$nCount_ATAC)
  equil_pt = equil_of_normals(my_mix$mu[mu_order[1]], my_mix$sigma[mu_order[1]], my_mix$lambda[mu_order[1]],
                              my_mix$mu[mu_order[2]], my_mix$sigma[mu_order[2]], my_mix$lambda[mu_order[2]])
  print(paste0("equil point is ", equil_pt) )

                      
  q11 = ggplot2::ggplot(observations) +
    ggplot2::geom_histogram(ggplot2::aes(x = nCount_ATAC), binwidth = 1, colour = "black", 
                            fill = "white") +
    ggplot2::stat_function(geom = "line", fun = function(x, mu, sigma, lam) { n*lam * dnorm(x, mu, sigma) },
                  args = list(my_mix$mu[mu_order[1]], my_mix$sigma[mu_order[1]], 
                              lam = my_mix$lambda[mu_order[1]]),
                  colour = "red", lwd = 1) +
    ggplot2::stat_function(geom = "line", fun = function(x, mu, sigma, lam) {n* lam * dnorm(x, mu, sigma) },
                  args = list(my_mix$mu[mu_order[2]], my_mix$sigma[mu_order[2]], 
                              lam = my_mix$lambda[mu_order[2]]),
                  colour = "blue", lwd = 1) +
    ggplot2::scale_x_continuous( limit = c(0, lower+100), oob = function(x, limits) x)+
    ggplot2::scale_y_continuous( limit = c(0, 3*n* my_mix$lambda[mu_order[2]]* dnorm(my_mix$mu[mu_order[2]], my_mix$mu[mu_order[2]], my_mix$sigma[mu_order[2]])), oob = function(x, limits) x)+
    ggplot2::ylab("Density") +
    ggplot2::xlab("Values") +
    ggplot2::ggtitle("Final GMM Fit")+
    ggplot2::geom_vline(xintercept = equil_pt,
                        # linetype="dotted",
                        color = "purple",
                        size=0.5)

  q2 = ggplot2::ggplot(observations, ggplot2::aes(x = nCount_ATAC)) +
    ggplot2::geom_histogram(binwidth =1) +
    mapply(
      function(mean, sd, lambda, n, binwidth) {
        ggplot2::stat_function(
          fun = function(x) {
            (stats::dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
          }
        )
      },
      mean = my_mix[["mu"]], #mean
      sd = my_mix[["sigma"]], #standard deviation
      lambda = my_mix[["lambda"]], #amplitude
      n = length(observations$nCount_ATAC), #sample size
      binwidth = 1  #binwidth used for histogram
    ) +
    ggplot2::scale_x_continuous( limit = c(0, lower+30), oob = function(x, limits) x)+
    ggplot2::scale_y_continuous( limit = c(0, 100000), oob = function(x, limits) x)+
    ggplot2::geom_vline(xintercept = lower,
               # linetype="dotted",
               color = "blue",
               size=0.5)+
    ggplot2::geom_vline(xintercept = barhop,
               # linetype="dotted",
               color = "red",
               size=0.5)


  q3 = ggplot2::ggplot(observations, ggplot2::aes(x = nCount_ATAC)) +
    ggplot2::geom_histogram(binwidth =1) +
    mapply(
      function(mean, sd, lambda, n, binwidth) {
        ggplot2::stat_function(
          fun = function(x) {
            (stats::dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
          }
        )
      },
      mean = my_mix[["mu"]], #mean
      sd = my_mix[["sigma"]], #standard deviation
      lambda = my_mix[["lambda"]], #amplitude
      n = length(observations$nCount_ATAC), #sample size
      binwidth = 1  #binwidth used for histogram
    ) +
    ggplot2::xlim(true_cells-my_mix$sigma[mu_order[3]], true_cells+my_mix$sigma[mu_order[3]])+
    ggplot2::ylim(0,75)

  # graphics::par(mfrow=c(1,3))
  # graphics::par(mar=c(5,5,4,1)+.1)
  # print(q1)
  # print(q2)
  # print(q3)

  egg::ggarrange(q1, q11, q2, q3 ,
            labels = c("A", "A2","B", "C"),
            ncol = 4, nrow = 1)

  return(c(lower, barhop))
}



# find point of equal probability

#' Find the point of equal probability between two Gaussian components of a Gaussian Mixture Model
#'
#' @param mu1 the mean of the first Gaussian distribution
#' @param sd1 the standard deviation of the first Gaussian distribution 
#' @param lam1 the mixing ratio of the first Gaussian distribution 
#' @param mu2 the mean of the second Gaussian distribution
#' @param sd2 the standard deviation of the second Gaussian distribution  
#' @param lam2 the mixing ratio of the second Gaussian distribution   
#'
#' @return
#'
equil_of_normals <- function(mu1, sd1, lam1, mu2, sd2, lam2){
  
  a = 1/sd1^2 - 1/sd2^2
  b = -2*( (mu1)/sd1^2 - (mu2)/sd2^2  )
  c = mu1^2/sd1^2 - mu2^2/sd2^2 + 2*log( lam2*sd1 / (sd2*lam1)  )
  discr = b^2 - 4*a*c
  
  if (a==0){
    return(-c/b)
  }
  
  if (discr<= 0 ){
    equil_pt = NULL
  } else {
    equil_pt = ( -b + base::sqrt(discr) ) / (2*a)
  }
  
  return(equil_pt)
  
}



# plot_counts <- function(nCount, max, min, vline1=NULL, vline2=NULL){
#   
#   
# }

