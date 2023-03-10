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








# sce <- Seurat::Read10X("/mnt/beegfs6/home3/ahringer/em613/repo_eD/eDv1/eD_multiome/data/input/valentina_FCA_GND10288180/raw_feature_bc_matrix")
#
# count_matrix_rna <- sce[["Gene Expression"]]
# count_matrix_atac <- sce[["Peaks"]]
#
# lower_rna = 200
# barhop_rna = 10
# lower_atac = 200
# barhop_atac = 7
#
# start_time <- Sys.time()
# eD.out_multi <- emptydrops_multiome(count_matrix_rna, lower_rna, barhop_rna, count_matrix_atac, lower_atac, barhop_atac )
# print("the number of cells detected is: ")
# print(sum(eD.out_multi$FDR_multi<0.001 & ! is.na(eD.out_multi$FDR_multi)))
# end_time <- Sys.time()
# print(end_time - start_time)
#

