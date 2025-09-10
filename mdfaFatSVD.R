mdfaFatSVD <- function(x, t) {
  n <- nrow(x)
  m <- ncol(x)
  p <- ncol(t)
  q <- p - m
  z <- x %*% t
  s <- svd(z)
  k <- s$u[, 1:m]
  l <- s$v[, 1:m]
  ll <- qr.Q(qr(cbind(l, matrix(rnorm(p * (p - m)), p, q))))
  kk <- qr.Q(qr(cbind(k, matrix(rnorm(n * (n - m)), n, n - m))))
  lperp <- ll[, -(1:m)]
  kperp <- kk[, -(1:m)]
  y <- tcrossprod(k, l) + tcrossprod(kperp[, 1:q], lperp)
  return(y)
}