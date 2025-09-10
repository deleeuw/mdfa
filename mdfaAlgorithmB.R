library(RSpectra)
source("mdfaAuxiliary.R")
source("mdfaProjections.R")

mdfaAlgorithmB <- function(cmat,
                           told,
                           proj = mdfaCFAProjection,
                           itmax = 1000,
                           eps = 1e-10,
                           verbose = TRUE) {
  itel <- 1
  m <- nrow(told)
  p <- ncol(told)
  q <- p - m
  ssqc <- sum(diag(cmat))
  ssqa <- sum(told^2)
  tct <- crossprod(told, cmat %*% told)
  ect <- eigs_sym(tct, m)
  fold <- ssqc + ssqa - 2 * sum(sqrt(ect$values))
  repeat {
    evc <- ect$vectors
    eva <- sqrt(ect$values)
    mva <- evc %*% diag(1 / eva) %*% t(evc)
    tnew <- proj(cmat %*% told %*% mva)
    tct <- crossprod(tnew, cmat %*% tnew)
    ect <- eigs_sym(tct, m)
    ssqa <- sum(tnew^2)
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

