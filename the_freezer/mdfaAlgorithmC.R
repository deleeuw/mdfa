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
  # yold <- projy(xmat %*% told, rank = m)
  # yold <- mdfaCompleteY(xmat, yold, told)
  yold <- matrix(rnorm(n * p), n, p)
  yold <- apply(yold, 2, function(x) x - mean(x))
  yold <- qr.Q(qr(yold))
  fold <- yold[, 1:q]
  uold <- yold[, -(1:q)]
  resi <- xmat - tcrossprod(fold, aold) - uold %*% dold
  sold <- sum(resi^2)
  itel <- 1
  repeat {
    dtil <- crossprod(uold, xmat)
    fnew <- projy((xmat - uold %*% dtil) %*% aold, rank = q)
    resi <- xmat - tcrossprod(fnew, aold) - uold %*% dold
    sfrf <- sum(resi^2)
    atil <- crossprod(xmat, fnew)
    anew <- aproj(atil)
    resi <- xmat - tcrossprod(fnew, anew) - uold %*% dold
    sfra <- sum(resi^2)
    mmat <- (xmat - tcrossprod(fnew, anew)) %*% dold
    unew <- projy((xmat - tcrossprod(fnew, anew)) %*% dold, rank = m)
    resi <- xmat - tcrossprod(fnew, anew) - unew %*% dold
    sfru <- sum(resi^2)
    dnew <- dproj(crossprod(unew, xmat))
    resi <- xmat - tcrossprod(fnew, anew) - unew %*% dnew
    sfrd <- sum(resi^2)
    if (verbose) {
      cat(
        "itel",
        formatC(itel, format = "d", width = 4),
        "sold",
        formatC(sold, digits = 10, format = "f"),
        "sfrf",
        formatC(sfrf, digits = 10, format = "f"),
        "sfra",
        formatC(sfra, digits = 10, format = "f"),
        "sfru",
        formatC(sfru, digits = 10, format = "f"),
        "sfrd",
        formatC(sfrd, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((sfrf - sfrd) < eps)) {
      break
    }
    itel <- itel + 1
    aold <- anew
    dold <- dnew
    uold <- unew
    fold <- fnew
    sold <- sfrd
  }
  return(list(
    ymat = cbind(unew, fnew),
    tmat = cbind(anew, dnew),
    loss = sfrd,
    itel = itel
  ))
}


