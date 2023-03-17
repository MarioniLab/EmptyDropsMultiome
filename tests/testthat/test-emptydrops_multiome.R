# test output dataframe is not empty when user specifies lowers and barhops and ATAC matrix has fewer droplets than RNA matrix
set.seed(0)
count_matrix <- DropletUtils:::simCounts()
colnames(count_matrix) <- paste0("C", as.character(seq(dim(count_matrix)[2]) ) )
count_matrix_rna <- as.matrix(count_matrix)
set.seed(1)
count_matrix1 <- DropletUtils:::simCounts()
colnames(count_matrix1) <- paste0("C", as.character(seq(dim(count_matrix1)[2]) ) )
count_matrix_atac <- as.matrix(count_matrix1)[, seq(1,11090)]

lower_rna = 200
barhop_rna = 10
lower_atac = 200
barhop_atac = 7

eD.out_multi <- emptydrops_multiome(count_matrix_rna, lower_rna, barhop_rna, count_matrix_atac, lower_atac, barhop_atac )
print("the number of cells detected is: ")
print(sum(eD.out_multi$FDR_multi<0.001 & ! is.na(eD.out_multi$FDR_multi)))

test_that("dataframe is not empty", {
  expect_equal(sum(eD.out_multi$Total_RNA) > 0, TRUE)
  expect_equal(sum(eD.out_multi$Total_chromatin) > 0, TRUE)
  expect_equal(sum(eD.out_multi$FDR_multi) > 0, TRUE)
})




# test output dataframe is not empty when lowers and barhops are inferred automatically and ATAC matrix has fewer droplets than RNA matrix
lower_rna = NULL
barhop_rna = NULL
lower_atac = NULL
barhop_atac = NULL

eD.out_multi <- emptydrops_multiome(count_matrix_rna, lower_rna, barhop_rna, count_matrix_atac, lower_atac, barhop_atac )
print("the number of cells detected is: ")
print(sum(eD.out_multi$FDR_multi<0.001 & ! is.na(eD.out_multi$FDR_multi)))

test_that("dataframe is not empty", {
  expect_equal(sum(eD.out_multi$Total_RNA) > 0, TRUE)
  expect_equal(sum(eD.out_multi$Total_chromatin) > 0, TRUE)
  expect_equal(sum(eD.out_multi$FDR_multi) > 0, TRUE)
})


# when there are only 2 droplets in each of the count matrices expect error 
tiny_count_matrix_rna <- count_matrix_rna[,1:2]
tiny_count_matrix_atac <- count_matrix_atac[,1:2]

lower_rna = NULL
barhop_rna = NULL
lower_atac = NULL
barhop_atac = NULL

test_that("when there are only two droplets expect error", {
  expect_error( emptydrops_multiome(tiny_count_matrix_rna, lower_rna, barhop_rna, tiny_count_matrix_atac, lower_atac, barhop_atac )  )
})




