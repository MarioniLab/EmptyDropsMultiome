
pvalues = c(1,1,1)

test_that("aggregating 1's gives 1", {
  expect_equal(fisher_v2(pvalues) ,1)
})

pvalues = c(0.3457261,0.187543,10e-319)

test_that("aggregating anything with 0 gives 0", {
  expect_equal(fisher_v2(pvalues), 0)
})


pvalues = c(0.3457261,0.187543, NA)

test_that("aggregating anything with NA gives NA", {
  expect_equal(fisher_v2(pvalues), NA)
})


#
# pvalues = c(1,1,1)
#
# test_that("aggregating 1's gives 1", {
#   expect_equal(ave_aggr(pvalues,2) ,1)
#   expect_equal(ave_aggr(pvalues,-1) ,1)
#   expect_equal(ave_aggr(pvalues,1) ,1)
# })
#


