library(RSpectra)
source("mdfaAuxiliary.R")
source("mdfaProjections.R")

mdfaAlgorithmW <- function(cmat,
                           wmat,
                           told,
                           projw = mdfaCFAWProjection,
                           itmax = 1000,
                           eps = 1e-10,
                           verbose = FALSE) {
  itel <- 1
  m <- nrow(told)
  p <- ncol(told)
  q <- p - m
  ssqc <- sum(wmat * cmat)
  ssqa <- sum(wmat * tcrossprod(told))
  tct <- crossprod(told, wmat %*% cmat %*% wmat %*% told)
  ect <- eigs_sym(tct, m)
  fold <- ssqc + ssqa - 2 * sum(sqrt(ect$values))
  repeat {
    evc <- ect$vectors
    eva <- sqrt(ect$values)
    mva <- evc %*% diag(1 / eva) %*% t(evc)
    tnew <- projw(cmat %*% wmat %*% told %*% mva, wmat)
    tct <- crossprod(tnew, wmat %*% cmat %*% wmat %*% tnew)
    ect <- eigs_sym(tct, m)
    ssqa <- sum(wmat * tcrossprod(tnew))
    fnew <- ssqc + ssqa - 2 * sum(sqrt(ect$values))
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
    fold <- fnew
  }
  return(list(
    tmat = tnew,
    loss = fnew,
    itel = itel
  ))
}

