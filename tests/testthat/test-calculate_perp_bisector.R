
point1 = c(0,1)
point2 = c(0,2)

line = perp_bisector(point1, point2)

test_that("horizontal lines are correct", {
  expect_equal(line, c(1.5, 0))
})

point1 = c(1,0)
point2 = c(2,0)

#line = perp_bisector(point1, point2)

# test_that("vertical lines lead to warning about default slope", {
#   expect_warning(perp_bisector(point1, point2))
# })

test_that("vertical lines are correct", {
  expect_equal(perp_bisector(point1, point2), c(1.5, 10^10))
})


