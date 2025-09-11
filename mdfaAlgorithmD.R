library(RSpectra)
source("mdfaAuxiliary.R")
source("mdfaProjections.R")

mdfaAlgorithmD <- function(x, yold, proj = mdfaCFAProjection, itmax = 1000, eps = 1e-10, verbose = TRUE) {
  told <- proj(crossprod(x, yold))
  zold <- x %*% told
  itel <- 1
  m <- ncol(x)
  resi <- x - tcrossprod(yold, told)
  fold <- sum(resi^2)
  repeat{
  ynew <- yold
  for (j in 1:m) {
    zj <- zold[, j]
    ynew[, j] <- 0
    yaux <- zj - ynew %*% drop(zj %*% ynew)
    ynew[, j] <- yaux / sqrt(sum(yaux ^ 2))
  }
  tnew <- proj(crossprod(x, ynew))
  znew <- x %*% tnew
  resi <- x - tcrossprod(ynew, tnew)
  fnew <- sum(resi^2)
  if (verbose) {
    cat(
      "itel",
      formatC(itel, format = "d", width = 4),
      "fold",
      formatC(fold, digits = 10, format = "f"),
      "fnew",
      formatC(fnew, digits = 10, format = "f"),
      "\n"
    )
  }
  if ((itel == itmax) || ((fold - fnew) < eps)) {
    break
  }
  itel <- itel + 1
  told <- tnew
  zold <- znew
  yold <- ynew
  fold <- fnew
  }
  return(list(
    y = ynew,
    t = tnew,
    loss = fnew,
    itel = itel
  ))
}