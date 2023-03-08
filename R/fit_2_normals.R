#' Title
#'
#' @param count_matrix matrix of counts
#'
#' @return
#' @export
#'
#' @examples
fit_2_normals <- function(count_matrix){

  # exp_counts_per_cell = unname(colSums(srat_029_all@assays[["RNA"]]@counts))
  #
  # my_mix <- normalmixEM(exp_counts_per_cell, mu=c(0,500,1000), epsilon=1e-09)
  # print(paste0("the mu's are : ", my_mix$mu))
  # my_mix$sigma
  #
  # mu_order = order(my_mix$mu)
  # mu_order
  #
  # lower = my_mix$mu[mu_order[2]]+2*my_mix$sigma[mu_order[2]]
  # lower
  # barhop = my_mix$mu[mu_order[1]]+4*my_mix$sigma[mu_order[1]]
  # barhop
  # observations = data.frame( "value" = exp_counts_per_cell )
  #
  # q1 = ggplot(observations, aes(x = value)) +
  #   geom_histogram(binwidth =1) +
  #   mapply(
  #     function(mean, sd, lambda, n, binwidth) {
  #       stat_function(
  #         fun = function(x) {
  #           (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
  #         }
  #       )
  #     },
  #     mean = my_mix[["mu"]], #mean
  #     sd = my_mix[["sigma"]], #standard deviation
  #     lambda = my_mix[["lambda"]], #amplitude
  #     n = length(observations$value), #sample size
  #     binwidth = 1  #binwidth used for histogram
  #   ) +
  #   xlim(0, 2000)+
  #   ylim(0,30000)+
  #   geom_vline(xintercept = lower,
  #              # linetype="dotted",
  #              color = "blue",
  #              size=0.5)+
  #   geom_vline(xintercept = barhop,
  #              # linetype="dotted",
  #              color = "red",
  #              size=0.5)
  # q1
  #
  # q2 = ggplot(observations, aes(x = value)) +
  #   geom_histogram(binwidth =1) +
  #   mapply(
  #     function(mean, sd, lambda, n, binwidth) {
  #       stat_function(
  #         fun = function(x) {
  #           (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
  #         }
  #       )
  #     },
  #     mean = my_mix[["mu"]], #mean
  #     sd = my_mix[["sigma"]], #standard deviation
  #     lambda = my_mix[["lambda"]], #amplitude
  #     n = length(observations$value), #sample size
  #     binwidth = 1  #binwidth used for histogram
  #   ) +
  #   xlim(0, lower)+
  #   ylim(0,300000)
  # q2
  #
  #
  # q3 = ggplot(observations, aes(x = value)) +
  #   geom_histogram(binwidth =1) +
  #   mapply(
  #     function(mean, sd, lambda, n, binwidth) {
  #       stat_function(
  #         fun = function(x) {
  #           (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
  #         }
  #       )
  #     },
  #     mean = my_mix[["mu"]], #mean
  #     sd = my_mix[["sigma"]], #standard deviation
  #     lambda = my_mix[["lambda"]], #amplitude
  #     n = length(observations$value), #sample size
  #     binwidth = 1  #binwidth used for histogram
  #   ) +
  #   xlim(1000, 3000)+
  #   ylim(0,100)
  #
  # print(q1 | q2 |q3 )
  # q1
  # q2
  # q3
  #
  # return(c(lower, barhop))
}
