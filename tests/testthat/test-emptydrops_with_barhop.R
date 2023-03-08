
lower = 400
barhop_end = 10
set.seed(0)
count_matrix <- DropletUtils:::simCounts()
colnames(count_matrix) <- paste0("C", as.character(seq(dim(count_matrix)[2]) ) )
#count_matrix <- as.matrix(count_matrix)
out = emptydrops_with_barhop(count_matrix, lower, barhop_end, niters=10000)

test_that("out has 6 attributes and 11100 rows", {
  expect_equal(length(out), 6)
  expect_equal(out@nrows, 11100)
})


test_that("the non NA FDRs are as many as the droplets with library size smaller than or equal to lower", {
  expect_equal(sum(is.na(out@listData[["FDR"]])), sum(Matrix::colSums(count_matrix)<=lower))
})

test_that("the number of droplets with barhop_amb_tent tags are as many as they should", {
  expect_equal(  sum(out@listData[["barhop_amb_tent"]]==1)  ,   sum(Matrix::colSums(count_matrix)<=lower & Matrix::colSums(count_matrix)>=barhop_end) )
  expect_equal(  sum(out@listData[["barhop_amb_tent"]]==2)  ,   sum( Matrix::colSums(count_matrix)<barhop_end) )
  expect_equal(  sum(out@listData[["barhop_amb_tent"]]==3)  ,   sum( Matrix::colSums(count_matrix)>lower) )
})


#
# # Identify likely cell-containing droplets.
# out <- emptyDrops(my.counts)
# out
#
# is.cell <- out$FDR <= 0.001
# sum(is.cell, na.rm=TRUE)
#
# # Subsetting the matrix to the cell-containing droplets.
# # (using 'which()' to handle NAs smoothly).
# cell.counts <- my.counts[,which(is.cell),drop=FALSE]
# dim(cell.counts)
#
# # Check if p-values are lower-bounded by 'niters'
# # (increase 'niters' if any Limited==TRUE and Sig==FALSE)
# table(Sig=is.cell, Limited=out$Limited)
#
#
# test_that("multiplication works", {
#   expect_equal(2 * 2, 4)
# })
