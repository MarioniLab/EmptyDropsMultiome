# simulate scRNA reads
hops = 500
soups = 100
cells = 100
mean_hops = 7
mean_soups = 300
mean_cells = 1300
atac_downsampling = 0.8
nCount = round(c(  rnorm(hops, mean=mean_hops, sd=1), rnorm(soups, mean=mean_soups, sd=50) , rnorm(cells, mean=mean_cells, sd=150) ))
out = fit_3_normals(nCount)

# simulate scATAC reads
hops = 500
soups = 100
cells = 100
mean_hops = 7
mean_soups = 300
mean_cells = 20000
nCount_atac = round(c(  rnorm(hops, mean=mean_hops, sd=1), rnorm(soups, mean=mean_soups, sd=50) , rnorm(cells, mean=mean_cells, sd=150) ))
out_atac = fit_3_normals_atac(nCount_atac)

test_that("RNA: lower is bigger than barhop_end", {
  expect_equal(out[1]>out[2], TRUE)
})
test_that("RNA: barhop_end is between the mean of barhops and the mean of ambient", {
  expect_equal(mean_hops<out[2] & mean_soups>out[2] , TRUE)
})
test_that("RNA: lower is larger than the mean of ambient", {
  expect_equal(mean_hops<out[2] & mean_soups>out[2] , TRUE)
})

test_that("lower atac is bigger than barhop_end", {
  expect_equal(out_atac[1]>out_atac[2], TRUE)
})

test_that("ATAC: barhop_end is between the mean of barhops and the mean of ambient", {
  expect_equal(mean_hops*atac_downsampling<out_atac[2] & mean_soups*atac_downsampling>out_atac[2] , TRUE)
})



# Cause bad convergence in ATAC and RNA
hops = 1
soups = 1
cells = 1
max_hops = 10
max_soups = 50
mean_cells = 20000
nCount_atac = round(c(  runif(hops, min=0, max=max_hops), runif(soups, min=max_hops, max=max_soups) , rnorm(cells, mean=mean_cells, sd=150) ))
out_atac = fit_3_normals_atac(nCount_atac)
nCount = round(c(  runif(hops, min=0, max=max_hops), runif(soups, min=max_hops, max=max_soups) , rnorm(cells, mean=mean_cells, sd=150) ))
out = fit_3_normals(nCount)
test_that("in case of bad convergence test default parameters are activated", {
  expect_equal( out_atac[1], 160 )
  expect_equal( out_atac[2], 4.625116 )
  expect_equal( out[1], 290.00000 )
  expect_equal( out[2], 29.567552 )
})



