
n=1000
m=200
lower = 100
barhop_end = 10
rna_count  = c( runif(n, 0, 500), runif(n, 700, 5000), runif(m, 0, 500), runif(m, 700, 5000) )
atac_count = c( runif(n, 0, 200), runif(n, 400, 6000), runif(m, 700, 5000), runif(m, 0, 500) )
df <- S4Vectors::DataFrame( "Total_RNA"= rna_count, "Total_chromatin"= atac_count,
                  "excluded"= rep(F,2*(n+m)),
                  "is_cell"= rep(0,2*(n+m)),
                  "FDR_multi"=runif(2*(n+m), 0, 0.002),
                  "barhop_amb_tent_RNA" =  as.numeric(rna_count<barhop_end) + as.numeric(rna_count<=lower) + as.numeric(rna_count>lower) *3,
                  "barhop_amb_tent_chromatin" =  as.numeric(atac_count<barhop_end) + as.numeric(atac_count<=lower) + as.numeric(atac_count>lower) *3
                  )

rownames(df)<- as.character(seq_along(rna_count))
#df_out <- call_cells(df)
df_out <- accept_k_means(df, lower, lower)



test_that("all droplets above ambiguous have been labelled as cells", {
  expect_equal(sum(df_out$FDR_multi[df_out$above_ambiguous==1]<0.001  ) , sum(df_out$above_ambiguous))
})

test_that("all droplets below barhop_end have been labelled as non cells", {
  expect_equal(sum(df_out$FDR_multi[df_out$Total_RNA<barhop_end]<0.001  ) , 0)
})

test_that("all droplets below barhop_end have been labelled as non cells", {
  expect_equal(sum(df_out$FDR_multi[df_out$Total_chromatin<barhop_end]<0.001  ) , 0)
})

