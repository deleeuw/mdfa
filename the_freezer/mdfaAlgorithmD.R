library(RSpectra)
source("mdfaAuxiliary.R")
source("mdfaProjections.R")

mdfaAlgorithmD <- function(xmat,
                           told,
                           proj = mdfaCFAProjection,
                           itmax = 1000,
                           eps = 1e-10,
                           verbose = TRUE) {
  itel <- 1
  m <- ncol(xmat)
  p <- ncol(told)
  q <- m - p
  n <- nrow(xmat)
  zold <- xmat %*% told
  yold <- projy(zold, rank = m)
  yold <- mdfaCompleteY(xmat, yold, told)
  resi <- xmat - tcrossprod(yold, told)
  fold <- sum(resi^2)
  repeat {
    ynew <- yold
    for (j in 1:p) {
      zolj <- zold[, j]
      ynew[, j] <- 0
      yaux <- zolj - ynew %*% drop(zolj %*% ynew)
      ynew[, j] <- yaux / sqrt(sum(yaux^2))
    }
    tnew <- proj(crossprod(xmat, ynew))
    znew <- xmat %*% tnew
    resi <- xmat - tcrossprod(ynew, tnew)
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
    ymat = ynew,
    tmat = tnew,
    loss = fnew,
    itel = itel
  ))
}