
hops = 500
soups = 100
cells = 100
mean_hops = 7
mean_soups = 300
mean_cells = 1300
atac_downsampling = 0.8

nCount = round(c(  rnorm(hops, mean=mean_hops, sd=1), rnorm(soups, mean=mean_soups, sd=50) , rnorm(cells, mean=mean_cells, sd=150) ))
out = fit_3_normals(nCount)


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



# # more overlapping Gaussian
# 
# hops = 1000
# soups = 100
# cells = 100
# mean_hops = 1
# mean_soups = 35
# mean_cells = 1300
# # hops = 5000
# # soups = 1000
# # cells = 1000
# atac_downsampling = 1
# 
# nCount = round(c(  rnorm(hops, mean=mean_hops, sd=10), rnorm(soups, mean=mean_soups, sd=30) , rnorm(cells, mean=mean_cells, sd=150) ))
# out_atac = fit_3_normals_atac(nCount)
# 
# equil = equil_of_normals(0.539105692069073, 10.8688383006814, 0.714113582232385, 31.6969572827669, 10.3549610529109, 0.142857142857143)
# equil
# 
# equil = equil_of_normals(100, 10, 0.5, 200, 10, 0.5)
# test_that("", {
#   expect_equal(equil, 150)
# })
# equil = equil_of_normals(100,10.1,0.5, 200, 10, 0.5)
# test_that("", {
#   expect_equal(equil, 150)
# })


