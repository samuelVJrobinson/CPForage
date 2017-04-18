context("Alpha function")

test_that("alpha is between 0 and 1",{
  expect_equal(alpha(0.05,1,14.35),0.01435)
  expect_equal(alpha(0.05,10,14.35),0.1435)
  expect_equal(alpha(0.05,50,14.35),0.7175)
})
