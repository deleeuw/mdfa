library(RSpectra)
source("mdfaAuxiliary.R")
source("mdfaProjections.R")

mdfaAlgorithmC <- function(xmat,
                           told,
                           aproj = proja,
                           dproj = projd,
                           itmax = 1000,
                           eps = 1e-10,
                           verbose = TRUE) {
  p <- ncol(told)
  m <- ncol(xmat)
  q <- p - m
  n <- nrow(xmat)
  aold <- told[, 1:q]
  dold <- told[, -(1:q)]
  yold <- xmat %*% told
  yold <- projy(xmat %*% told, rank = m)
  fold <- qr.Q(qr(yold[, 1:q]))
  uold <- qr.Q(qr(yold[, -(1:q)]))
  resi <- xmat - tcrossprod(fold, aold) - uold %*% dold
  sold <- sum(resi^2)
  itel <- 1
  repeat {
    dtil <- crossprod(uold, xmat)
    fnew <- projy((xmat - uold %*% dtil) %*% aold, rank = q)
    atil <- crossprod(xmat, fnew)
    anew <- aproj(atil)
    resi <- xmat - tcrossprod(fnew, anew) - uold %*% dold
    smid <- sum(resi^2)
    unew <- projy((xmat - tcrossprod(fnew, anew)) %*% dold, rank = m)
    dtil <- crossprod(unew, xmat)
    dnew <- dproj(dtil)
    resi <- xmat - tcrossprod(fnew, anew) - unew %*% dnew
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
    ymat <- cbind(unew, fnew),
    tmat <- cbind(anew, dnew),
    loss = snew,
    itel = itel
  ))
}


