
e.out1 <- data.frame( "PValue"=c(1,1,1,0,0,0,0.5,0.5)   )
rownames(e.out1) <- as.character(c(seq_along(e.out1$PValue)) )

e.out2 <- data.frame( "PValue"=c(1,1,0,0,0,1,0.5,   0.5,0.5,0.5)   )
rownames(e.out2) <- c(rownames(e.out1)[seq(1,7)], as.character(c(9,10,11)) )

e_multi.out <- aggr_p_values( e.out1, e.out2)

test_that("merging doesn't lose droplets", {
  expect_equal( sort(union(rownames(e.out1), rownames(e.out2) )) , rownames(e_multi.out)   )
})


# test_that("fischer+mean works", {
#   expect_equal(e_multi.out[as.character(c( seq(1,6), seq(8,11) )),]$PValue_multi ,  c(1,1,0,0,0,0,    NA,NA,NA,NA)   )
# })
