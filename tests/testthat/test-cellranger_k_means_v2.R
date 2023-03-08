
data <- data.frame( "dup" = 1 - c(T,T,F,F,T,F,F,F,T,F,T,T,F),
                    "is_cell" = c(1,1,0,0,1,0,0,0,0,0,0,1,0) )
corrected_is_cell = c(rep(1,8), c(0,0,0), c(1,1))

reconstructed_is_cell = transfer_labels_to_dedupl_v3(data)$is_cell

test_that("label transfer works", {
  expect_equal(corrected_is_cell, reconstructed_is_cell)
})

n=1000
rna_count = c( runif(n, 0, 500), runif(n, 700, 5000) )
atac_count = c( runif(n, 0, 200), runif(n, 400, 6000) )
df <- data.frame( "rna_count"= rna_count, "atac_count"= atac_count,
              "excluded"= rep(F,2*n),
              "is_cell"= rep(0,2*n)    )
df_out <- call_cells(df)

test_that("we have 1000 cells", {
  expect_equal(sum(df_out[[1]]$is_cell), n)
})

test_that("we have 1000 empty droplets", {
  expect_equal(sum(! df_out[[1]]$is_cell), n)
})

