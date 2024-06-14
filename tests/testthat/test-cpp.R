test_that("test-cpp", {
  expect_equal(rowSumsC(matrix(c(1,2,3, 11,12,13), nrow = 2, ncol = 3, byrow = TRUE)), c(6, 36))
})
