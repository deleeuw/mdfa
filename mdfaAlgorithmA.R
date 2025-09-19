library(RSpectra)
source("mdfaAuxiliary.R")
source("mdfaProjections.R")

mdfaAlgorithmA <- function(xmat,
                           told,
                           proj = mdfaCFAProjection,
                           itmax = 1000,
                           eps = 1e-10,
                           ortho = TRUE,
                           verbose = TRUE) {
  m <- nrow(told)
  yold <- projy(xmat %*% told, rank = m)
  resi <- xmat - tcrossprod(yold, told)
  fold <- sum(resi^2)
  itel <- 1
  repeat {
    tnew <- proj(crossprod(xmat, yold))
    ynew <- projy(xmat %*% tnew, rank = m)
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
    yold <- ynew
    fold <- fnew
  }
  if (ortho) {
    ynew <- mdfaCompleteY(xmat, ynew, tnew) 
  }
  return(list(
    ymat = ynew,
    tmat = tnew,
    loss = fnew,
    itel = itel
  ))
}
