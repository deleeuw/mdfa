
mdfaAlgorithmG <- function(cmat, told) {
  m <- nrow(told)
  p <- ncol(told)
  q <- p - m
  h <- optim(
    par = as.vector(told),
    fn = theFunc,
    gr = theGrad,
    method = "BFGS",
    control = list(trace = 6),
    cmat = cmat,
    nrow = m,
    ncol = p
  )
  tmat <- matrix(h$par, m, p)
  return(
    list(
      tmat = tmat,
      loss = h$value,
      counts = h$counts
    )
  )
}

theFunc <- function(told, cmat, nrow, ncol) {
  tcp <- sum(diag(cmat))
  tap <- sum(told^2)
  amt <- matrix(told, nrow, ncol)
  aca <- crossprod(amt, cmat %*% amt)
  func <- tcp + tap - 2 * sum(sqrt(svd(aca)$d[1:nrow]))
  return(func)
}

theGrad <- function(a, cmat, nrow, ncol) {
  amt <- matrix(a, nrow, ncol)
  aca <- crossprod(amt, cmat %*% amt)
  eca <- eigen(aca)
  evc <- eca$vectors[, 1:nrow]
  eva <- diag(1 / sqrt(eca$values[1:nrow]))
  mva <- evc %*% eva %*% t(evc)
  return(as.vector(a - cmat %*% amt %*% mva))
}




