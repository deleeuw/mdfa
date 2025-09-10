library(RSpectra)
source("mdfaAuxiliary.R")
source("mdfaProjections.R")

proja <- function(x) {
  return(x)
}

projd <- function(x) {
  return(diag(diag(x)))
}

projy <- function(x) {
  s <- svd(x)
  return(tcrossprod(s$u, s$v))
}

mdfaAlgorithmC <- function(x,
                           aold,
                           dold,
                           aproj = proja,
                           dproj = projd,
                           itmax = 1000,
                           eps = 1e-10,
                           verbose = TRUE) {
  q <- ncol(aold)
  m <- nrow(aold)
  p <- q + m
  told <- cbind(aold, dold)
  yold <- x %*% told
  yold <- projy(x %*% told)
  fold <- qr.Q(qr(yold[, 1:q]))
  uold <- qr.Q(qr(yold[, -(1:q)]))
  resi <- x - tcrossprod(fold, aold) - uold %*% dold
  sold <- sum(resi^2)
  itel <- 1
  repeat {
    dtil <- crossprod(uold, x)
    fnew <- projy((x - uold %*% dtil) %*% aold)
    atil <- crossprod(x, fnew)
    anew <- aproj(atil)
    resi <- x - tcrossprod(fnew, anew) - uold %*% dold
    smid <- sum(resi^2)
    unew <- projy((x - tcrossprod(fnew, anew)) %*% dold)
    dtil <- crossprod(x, unew)
    dnew <- dproj(dtil)
    resi <- x - tcrossprod(fnew, anew) - unew %*% dnew
    snew <- sum(resi^2)
    if (verbose) {
      cat(
        "itel",
        formatC(itel, format = "d", width = 4),
        "sold",
        formatC(sold, digits = 10, format = "f"),
        "smid",
        formatC(smid, digits = 10, format = "f"),
        "snew",
        formatC(snew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    itel <- itel + 1
    aold <- anew
    dold <- dnew
    uold <- unew
    fold <- fnew
    sold <- snew
  }
  return(list(
    a = anew,
    d = dnew,
    u = unew,
    f = fnew,
    loss = snew,
    itel = itel
  ))
}


