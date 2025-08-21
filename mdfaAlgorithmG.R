
mdfaAlgorithmG <- function(c, a, t) {
  nrow <- nrow(a)
  ncol <- ncol(a)
  ncom <- ncol - nrow
  h <- optim(
    par = as.vector(a),
    fn = theFunc,
    gr = theGrad,
    method = "BFGS",
    control = list(trace = 6),
    cp = c,
    tp = t,
    nrow = nrow,
    ncol = ncol
  )
  hmat <- matrix(h$par, nrow, ncol)
  loadings <- hmat[, 1:ncom]
  uniquenesses <- diag(hmat[, -(1:ncom)]^2)
  mat <- crossprod(loadings, diag(1 / uniquenesses) %*% loadings)
  mvc <- eigen(mat)$vectors
  loadings <- loadings %*% mvc
  return(
    list(
      loadings = loadings,
      uniquenesses = uniquenesses,
      loss = h$value,
      counts = h$counts
    )
  )
}

theFunc <- function(a, cp, tp, nrow, ncol) {
  tcp <- sum(diag(cp))
  tap <- sum(a^2)
  amt <- matrix(a, nrow, ncol)
  aca <- crossprod(amt, cp %*% amt)
  func <- tcp + tap - 2 * sum(sqrt(svd(aca)$d[1:nrow]))
  return(func)
}

theGrad <- function(a, cp, tp, nrow, ncol) {
  amt <- matrix(a, nrow, ncol)
  aca <- crossprod(amt, cp %*% amt)
  eca <- eigen(aca)
  evc <- eca$vectors[, 1:nrow]
  eva <- diag(1 / sqrt(eca$values[1:nrow]))
  mva <- evc %*% eva %*% t(evc)
  return(as.vector(tp * (a - cp %*% amt %*% mva)))
}




