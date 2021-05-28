test_that("spatTree builds a tree", {
  y <- c(rep(0,5),rep(1,5))
  x <- cbind(c(1:10),rep(c(1,2),5))
  y1 <- spatTree(y,x,i.Sig=diag(10),m=2)
  expect_equal(c(1,5.5,1,2.5),as.numeric(unlist(y1))[1:4])
})
