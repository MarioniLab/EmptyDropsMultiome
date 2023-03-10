#' Finds the droplets that are in 1) the top 1% of RNA counts and in the top 0.1% of atac counts, and then
#' calculates their mean RNA counts and mean ATAC counts.
#'
#' @param data Dataframe with data$rna being the number of RNA counts, and data$atac being the number of ATAC counts.
#'
#' @return A real vector with the two mean values.
#'
find_tip <- function(data){
  # data: dataframe with data$rna being the number of RNA counts, and data$atac being the number of ATAC counts
  # output: center coordinates of the center of the droplets in the top percentile
  p99_rna = stats::quantile(data$Total_RNA, probs = c(0.99), na.rm=TRUE )
  p999_atac = stats::quantile(data$Total_chromatin, probs = c(0.999), na.rm=TRUE )
  print("the tip has rna count > ")
  print(p99_rna)
  print("the tip has atac count > ")
  print(p999_atac)
  tip = data$Total_RNA > p99_rna & data$Total_chromatin > p999_atac

  if (sum(tip)==0){
    p95_rna = stats::quantile(data$Total_RNA, probs = c(0.95), na.rm=TRUE )
    p95_atac = stats::quantile(data$Total_chromatin, probs = c(0.95), na.rm=TRUE )
    tip = data$Total_RNA > p95_rna & data$Total_chromatin > p95_atac

  }

  tip_center <- c( mean(data$Total_RNA[tip]),   mean(data$Total_chromatin[tip]) )

  return(tip_center)

}

#' Calculates the distance between a point and a line.
#'
#' @param equation A vector c(intercept, slope) that defines a line in 2d Euclidean space.
#' @param point A 2d vector of coordinates.
#'
#' @return A single real number.
#'
dist_point_to_line <- function(equation, point){
  # equation: numeric vector for intercept=equation[1] and slope=equation[2]
  # point: numeric vector for coordinates of point
  # output: dist, the distance of the point to the line (in 2d)

  intercept = equation[1]
  slope = equation[2]

  # if line has the form a*x+b*y+c=0
  a = slope
  b = -1
  c = intercept

  # point coordinates
  x0 = point[1]
  y0 = point[2]

  dist = abs(a*x0+b*y0+c)/sqrt(a^2+b^2)
  return(dist)

}

#' Calculate the equation of a line that is parallel to a given line and a given distance from it.
#'
#' @param equation A 2d vector c(intercept, slope) defining a line in 2d Euclidean space.
#' @param dist A real number specifying the distance between the given line and the line to be calculated.
#'
#' @return A 2d vector c(intercept, slope) defining a line in 2d Euclidean space
#'
calc_intercept_of_parallel_line <- function(equation, dist){
  # equation: numeric vector for intercept=equation[1] and slope=equation[2]
  # dist: real number for the distance between the given line and the line to be calculated
  # output: equation2, equation for the line to be calculated

  intercept1 = equation[1]
  slope = equation[2]

  intercept2 = intercept1 + dist / sqrt(slope^2 + 1)

  return( c( intercept2, slope) )

}

#' Calculate a line that is parallel to the k-means line and a distance from it that is equal to 66% of the distance
#' between the droplets with the biggest library size and the k-means line. The area bound by the k-means line and this
#' line is referred to as ambiguous area and in it droplets are rejected by thresholding FDR_multi.
#' Above this ambiguous area, droplets are assumed to contain nuclei/cells and their FDR_multi is set to 0 by default.
#'
#' @param data Dataframe of droplets with columns Total_RNA and Total_chromatin.
#' @param equation 2d vector c(intercept, slope) defining the k-means line.
#'
#' @return 2d vector c(intercept, slope) defining the k-means line.
#'
calc_ambiguous_above <- function(data, equation){
  # given the data, calculate the ambiguous above area (specifically the line outlining its upper bound) (i.e. the area just north of the cellranger-arc line)
  tip_center <- find_tip(data)

  print("tip_center")
  print(tip_center)

  # log transform the coordinates
  tip_center = c( log10(tip_center[1]+0.1), log10(tip_center[2]+0.1) )

  dist = dist_point_to_line(equation, tip_center) / 1.5

  equation_parallel = calc_intercept_of_parallel_line(equation, dist)

  return(equation_parallel)


}




#' Calculate a line that is parallel to a given line and passes through the point (x, y).
#'
#' @param x A real number for the x coordinate of the point.
#' @param y A real number for the y coordinate of the point.
#' @param equation A 2d vector c(intercept, slope) defining a line in 2d Euclidean space.
#'
#' @return A 2d vector c(intercept, slope) defining a line in 2d Euclidean space.
#'
calc_lower_line <- function(x, y, equation){
  # calculate a line equation_lower which is parallel to equation and goes through the point (x, y)

  slope = equation[2]
  intercept = y - slope * x

  return( c( intercept, slope) )
}


