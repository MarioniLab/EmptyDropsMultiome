#' Fits a mixture of three normal distributions to histogram of RNA counts, to infer the parameters lower and barhop_end
#' which are used in emptydrops.
#'
#' @param exp_counts_per_cell a vector with RNA library size per barcode.
#' @param verbose if TRUE various intermediate steps and plots are printed.
#'
#' @return A 2d vector with the values lower and barhop_end for RNA.
#'
#' @export
#'
fit_3_normals <- function(exp_counts_per_cell, verbose=TRUE){

  my_mix <- tryCatch(
    {
      mixtools::normalmixEM(exp_counts_per_cell, mu=c(5,500,1200), sigma=c(5,50,800) , epsilon=1e-09)
    },
    error=function(cond) {
      message("Use default settings for barhop_end and lower")
      # Choose a default values in case of convergence error
      return(data.frame("mu"=c(25, 200, 1000  ), "sigma"=c(1, 60, 300  ), "lambda"= c(1/1.15, 0.10/1.15, 0.05/1.15)))
    },
    warning=function(cond) {},
    finally={}
  )    
  
  
  n = length(exp_counts_per_cell)
  mu_order = order(my_mix$mu)
  mu_barhop = my_mix$mu[mu_order[1]]
  mu_ambient = my_mix$mu[mu_order[2]]
  mu_true_cells = my_mix$mu[mu_order[3]]
  
  if (mu_barhop > 100 | mu_ambient >700 | mu_barhop>mu_ambient-15){
    my_mix <- data.frame("mu"=c(25, 200, 1000  ), "sigma"=c(1, 60, 300  ), "lambda"= c(1/1.15, 0.10/1.15, 0.05/1.15))
    n = length(exp_counts_per_cell)
    mu_order = order(my_mix$mu)
    mu_barhop = my_mix$mu[mu_order[1]]
    mu_ambient = my_mix$mu[mu_order[2]]
    mu_true_cells = my_mix$mu[mu_order[3]]
  }

  lower = my_mix$mu[mu_order[2]]+1.5*my_mix$sigma[mu_order[2]]
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

    plot_counts(exp_counts_per_cell, mu_barhop, my_mix$sigma[mu_order[1]], my_mix$lambda[mu_order[1]],
                mu_ambient, my_mix$sigma[mu_order[2]], my_mix$lambda[mu_order[2]], 
                vline1 =  equil_pt, vline2 = lower)
    
  }

  return(c(lower, equil_pt))
}







#' Fits a mixture of three normal distributions to histogram of ATAC counts, to infer the parameters lower and barhop_end
#' which are used in emptydrops.
#'
#' @param exp_counts_per_cell a vector with the ATAC library size per barcode.
#' @param verbose if TRUE various intermediate steps and plots are printed.
#'
#' @return A 2d vector with the values lower and barhop_end for ATAC.
#'
#' @export
#'
fit_3_normals_atac <- function(exp_counts_per_cell, verbose=TRUE){
  
  # t <- try(my_mix <- mixtools::normalmixEM(exp_counts_per_cell, lambda=c(6*10^5/(7*10^5), 8*10^4/(7*10^5), 10^4/(7*10^5) ), 
  #                                          mu=c(5,100,10000), sigma=c(5,50,800) , epsilon=1e-09))
  # print(inherits(t,"try-error"))
  # if( inherits(t,"try-error")  ){
  #   my_mix <- data.frame("mu"=c(1, 40, 800  ), "sigma"=c(1, 60, 600  ), "lambda"= c(1,0.10,0.05))
  # }
  # 
  
  my_mix <- tryCatch(
    {
      mixtools::normalmixEM(exp_counts_per_cell, lambda=c(6*10^5/(7*10^5), 8*10^4/(7*10^5), 10^4/(7*10^5) ), 
                            mu=c(5,100,10000), sigma=c(5,50,800) , epsilon=1e-09)
    },
    error=function(cond) {
      message("Use default settings for barhop_end and lower")
      # Choose a default values in case of convergence error
      return(data.frame("mu"=c(1, 40, 800  ), "sigma"=c(1, 60, 600  ), "lambda"= c(1,0.10,0.05)))
    },
    warning=function(cond) {},
    finally={}
  )    
  
  
  n = length(exp_counts_per_cell)
  mu_order = order(my_mix$mu)
  mu_barhop = my_mix$mu[mu_order[1]]
  mu_ambient = my_mix$mu[mu_order[2]]
  mu_true_cells = my_mix$mu[mu_order[3]]
  
  if (mu_barhop > 60 | mu_ambient >600 | mu_barhop>mu_ambient-15){
    my_mix <- data.frame("mu"=c(1, 40, 800  ), "sigma"=c(1, 60, 600  ), "lambda"= c(1,0.10,0.05))
    mu_order = order(my_mix$mu)
    mu_barhop = my_mix$mu[mu_order[1]]
    mu_ambient = my_mix$mu[mu_order[2]]
    mu_true_cells = my_mix$mu[mu_order[3]]
  }
  
  equil_pt = equil_of_normals(mu_barhop, my_mix$sigma[mu_order[1]], my_mix$lambda[mu_order[1]],
                              mu_ambient, my_mix$sigma[mu_order[2]], my_mix$lambda[mu_order[2]])
  lower = my_mix$mu[mu_order[2]]+2*my_mix$sigma[mu_order[2]]
  true_cells = my_mix$mu[mu_order[3]]
  
  
  if (verbose){
    
    print(paste0("the mu's are : ", my_mix$mu[1],", ", my_mix$mu[2],", ", my_mix$mu[3] ))
    my_mix$mu
    print(paste0("the sigma's are : ", my_mix$sigma[1],", ", my_mix$sigma[2],", ", my_mix$sigma[3] ))
    my_mix$sigma
    print(paste0("the lambdas are : ", my_mix$lambda[1],", ", my_mix$lambda[2],", ", my_mix$lambda[3] ))
    my_mix$lambda
    print(paste0("equil point is ", equil_pt) )
    print(paste0("the lower is : ", lower ))
    
    plot_counts(exp_counts_per_cell, mu_barhop, my_mix$sigma[mu_order[1]], my_mix$lambda[mu_order[1]],
                mu_ambient, my_mix$sigma[mu_order[2]], my_mix$lambda[mu_order[2]], 
                vline1 =  equil_pt, vline2 = lower)
    
  }
  

  return(c(lower, equil_pt))
}



# Find point of equal probability

#' Find the point of equal probability between two Gaussian components of a Gaussian Mixture Model
#'
#' @param mu1 the mean of the first Gaussian distribution
#' @param sd1 the standard deviation of the first Gaussian distribution 
#' @param lam1 the mixing ratio of the first Gaussian distribution 
#' @param mu2 the mean of the second Gaussian distribution
#' @param sd2 the standard deviation of the second Gaussian distribution  
#' @param lam2 the mixing ratio of the second Gaussian distribution   
#'
#' @return Real number for the point of intersection of the Gaussians
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



#' Plot histogram of counts for lower and barhop_end inference
#'
#' @param nCount vector with the RNA or ATAC counts
#' @param mu1 mu of first Gaussian
#' @param sigma1 sigma of first Gaussian
#' @param lam1 lambda of first Gaussian
#' @param mu2 mu of first Gaussian
#' @param sigma2 sigma of first Gaussian
#' @param lam2 lambda of first Gaussian
#' @param max upper limit of RNA range
#' @param vline1 x value of first vertical line
#' @param vline2 x value of second vertical line
#' @param label label for x axis
#'
#' @return None
#' @export
#'
plot_counts <- function(nCount, mu1, sigma1, lam1, mu2, sigma2, lam2, max=50000, vline1, vline2, label){
  
  observations = data.frame( "nCount" = nCount )
  
  n = length(observations$nCount)
  
  overview = ggplot2::ggplot(observations) + ggplot2::theme_bw() +
    ggplot2::geom_histogram(ggplot2::aes(x = nCount), binwidth = 1, colour = "black", 
                            fill = "black") +
    ggplot2::stat_function(geom = "line", fun = function(x, mu, sigma, lam) { n*lam * stats::dnorm(x, mu, sigma) },
                           args = list(mu1, sigma1, lam = lam1),
                           colour = "red", lwd = 1) +
    ggplot2::stat_function(geom = "line", fun = function(x, mu, sigma, lam) {n* lam * stats::dnorm(x, mu, sigma) },
                           args = list(mu2, sigma2, lam = lam2),
                           colour = "blue", lwd = 1) +
    ggplot2::scale_x_continuous( limit = c(0, vline2+100), oob = function(x, limits) x)+
    ggplot2::scale_y_continuous( limit = c(0, max), oob = function(x, limits) x)+
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("nCount") +
    ggplot2::ggtitle("Final GMM Fit")+
    ggplot2::geom_vline(xintercept = vline1,
                        # linetype="dotted",
                        color = "purple",
                        size=0.5)+
  ggplot2::geom_vline(xintercept = vline2,
                      # linetype="dotted",
                      color = "blue",
                      size=0.5)
  
  overview_nofit = ggplot2::ggplot(observations) + ggplot2::theme_bw() +
    ggplot2::geom_histogram(ggplot2::aes(x = nCount), binwidth = 1, colour = "black", 
                            fill = "black") +
    ggplot2::scale_x_continuous( limit = c(0, vline2+100), oob = function(x, limits) x)+
    ggplot2::scale_y_continuous( limit = c(0, max), oob = function(x, limits) x)+
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("nCount") +
    ggplot2::ggtitle("Final GMM Fit")+
    ggplot2::geom_vline(xintercept = vline1,
                        # linetype="dotted",
                        color = "purple",
                        size=0.5)+
    ggplot2::geom_vline(xintercept = vline2,
                        # linetype="dotted",
                        color = "blue",
                        size=0.5)

  
  zoomin = ggplot2::ggplot(observations) + ggplot2::theme_bw() + 
    ggplot2::geom_histogram(ggplot2::aes(x = nCount), binwidth = 1, colour = "black", 
                            fill = "black") +
    ggplot2::stat_function(geom = "line", fun = function(x, mu, sigma, lam) { n*lam * stats::dnorm(x, mu, sigma) },
                           args = list(mu1, sigma1, lam = lam1),
                           colour = "red", lwd = 1) +
    ggplot2::stat_function(geom = "line", fun = function(x, mu, sigma, lam) {n* lam * stats::dnorm(x, mu, sigma) },
                           args = list(mu2, sigma2, lam = lam2),
                           colour = "blue", lwd = 1) +
    ggplot2::scale_x_continuous( limit = c(0, vline2+100), oob = function(x, limits) x)+
    ggplot2::scale_y_continuous( limit = c(0, 3*n* lam2* stats::dnorm(mu2, mu2, sigma2)), oob = function(x, limits) x)+
    ggplot2::ylab("Frequency") +
    ggplot2::xlab("nCount") +
    ggplot2::ggtitle("Final GMM Fit")+
    ggplot2::geom_vline(xintercept = vline1,
                        # linetype="dotted",
                        color = "purple",
                        size=0.5)+
    ggplot2::geom_vline(xintercept = vline2,
                        # linetype="dotted",
                        color = "blue",
                        size=0.5)
  
  print(overview) 
  print(overview_nofit) 
  print(zoomin)
  
}

