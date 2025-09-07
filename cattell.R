data(cattell, package = "psych")

cattellProjection <- function(tmat) {
  tmat[-(1:2), 2] <- 0
  tmat[-(3:4), 3] <- 0
  tmat[-(5:6), 4] <- 0
  tmat[-(7:8), 5] <- 0
  tmat[-(9:12), 6] <- 0
  tmat[, 7:18] <- diag(diag(tmat[, 7:18]))
  return(tmat)
}

cattellSVDProjection <- function(tmat) {
  tmat[1:2, 2:6] <- rankOne(tmat[1:2, 2:6])
  tmat[3:4, 2:6] <- rankOne(tmat[3:4, 2:6])
  tmat[5:6, 2:6] <- rankOne(tmat[5:6, 2:6])
  tmat[7:8, 2:6] <- rankOne(tmat[7:8, 2:6])
  tmat[9:12, 2:6] <- rankOne(tmat[9:12, 2:6])
  tmat[, 7:18] <- diag(diag(tmat[, 7:18]))
  return(tmat)
}

cattellUnravel <- function(tmat) {
  s2 <- svd(tmat[1:2, 2:6])
  s3 <- svd(tmat[3:4, 2:6])
  s4 <- svd(tmat[5:6, 2:6])
  s5 <- svd(tmat[7:8, 2:6])
  s6 <- svd(tmat[9:12, 2:6])
  a <- matrix(0, 12, 6)
  a[, 1]<- tmat[, 1]
  r <- diag(6)
  d <- diag(tmat[, 7:18] ^ 2)
  a[1:2, 2] <- s2$u[, 1]
  a[3:4, 3] <- s3$u[, 1]
  a[5:6, 4] <- s4$u[, 1]
  a[7:8, 5] <- s5$u[, 1]
  a[9:12, 6] <- s6$u[, 1]
  r[2, -1] <- s2$d[1] * s2$v[, 1]
  r[3, -1] <- s3$d[1] * s3$v[, 1]
  r[4, -1] <- s4$d[1] * s4$v[, 1]
  r[5, -1] <- s5$d[1] * s5$v[, 1]
  r[6, -1] <- s6$d[1] * s6$v[, 1]
  s <- 1 / sqrt(rowSums(r ^ 2))
  r <- diag(s) %*% r
  a <- a %*% diag(1 / s)
  return(list(loadings = a, correlations = tcrossprod(r), uniquenesses = d))
}

groups <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5, 5)

cottellProjection <- function(tmat) {
  m <- nrow(tmat)
  p <- ncol(tmat)
  q <- p - m
  i <- 1:m
  templ <- ifelse(outer(groups, groups, "=="), 1, 0)
  templ <- templ * ifelse(outer(1:m, 1:m, ">="), 1, 0)
  tmat[, -(1:q)] <- templ * tmat[, -(1:q)]
  return(tmat)
}

rankOne <- function(x) {
  s <- svd(x)
  return(s$d[1] * outer(s$u[, 1], s$v[, 1]))
}

aeig <- eigen(cattell)
avec <- aeig$vectors[, 1:6]
aval <- diag(sqrt(aeig$values[1:6]))
acom <- avec %*% aval
catold <- cbind(acom, diag(sqrt(1 - rowSums(acom^2))))
ctmat <- cattellProjection(catold)
csmat <- cattellSVDProjection(catold)
avec <- aeig$vectors[, 1:2]
aval <- diag(sqrt(aeig$values[1:2]))
acom <- avec %*% aval
cemat <- cottellProjection(cbind(acom, matrix(rnorm(144), 12, 12)))

catemp <- ifelse(ctmat != 0, 1, 0)