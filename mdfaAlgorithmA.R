library(RSpectra)
source("mdfaAuxiliary.R")
source("mdfaProjections.R")

mdfaAlgorithmA <- function(x,
                           told,
                           proj = mdfaCFAProjection,
                           itmax = 1000,
                           eps = 1e-10,
                           ortho = FALSE,
                           verbose = TRUE) {
  m <- nrow(told)
  svxa <- svd(x %*% told, nu = m, nv = m)
  yold <- tcrossprod(svxa$u, svxa$v)
  resi <- x - tcrossprod(yold, told)
  fold <- sum(resi^2)
  itel <- 1
  repeat {
    tnew <- proj(crossprod(x, yold))
    svxa <- svd(x %*% tnew, nu = m, nv = m)
    ynew <- tcrossprod(svxa$u, svxa$v)
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
    yold <- ynew
    fold <- fnew
  }
  if (ortho) {
    kperp <- leftNullSpace(x %*% tnew)
    lperp <- leftNullSpace(t(x %*% tnew))
    ynew <- ynew + tcrossprod(kperp[, 1:q], lperp)
  }
  return(list(
    y = ynew,
    t = tnew,
    loss = fnew,
    itel = itel
  ))
}
