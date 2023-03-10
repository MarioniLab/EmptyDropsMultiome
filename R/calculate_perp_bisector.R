#' log-transform the data
#'
#' @param data Dataframe with column x and column y
#'
#' @return Dataframe with with log-transformed column values
#'
#'
log_transform_v2 <- function(data){
  # data: dataframe with column "atac_count" and "rna_count"
  # output: dataframe with log10 transformed columns
  data["x"] <- log10(data["x"]+0.1)
  data["y"] <- log10(data["y"]+0.1)
  return(data)
}

#' Calculate the perpendicular bisector line to the line segment defined by two points in 2d Euclidean space.
#'
#' @param center1 A point in 2d Euclidean space.
#' @param center2 A point in 2d Euclidean space.
#'
#' @return Two dimensional vector c(intercept, slope)
#'
perp_bisector <- function(center1, center2){
  # center1: numeric vector with position coordinates
  # center2: numeric vector with position coordinates
  # output: calculate slope (b) and intercept (a) of perpendicular bisector: y = b * x +a

  dy = center2[2]-center1[2]
  y_mean = (center2[2]+center1[2])/2

  dx = center2[1]-center1[1]
  x_mean = (center2[1]+center1[1])/2

  a = y_mean+ x_mean *(dx / dy)
  b = - dx / dy

  if (b==Inf || b==-Inf){
    #warning('slope became Infinite, default to slope=10^10')
    b=10^10
    a=x_mean
  }

  return( c(a,b) )


}




