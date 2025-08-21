mdfaAlgorithmA <- function(x,
                           aold,
                           template,
                           itmax = 1000,
                           eps = 1e-10,
                           verbose = TRUE) {
  m <- ncol(x)
  svxa <- svd(x %*% aold)
  yold <- tcrossprod(svxa$u[, 1:m], svxa$v[, 1:m])
  resi <- x - tcrossprod(yold, aold)
  fold <- sum(resi^2)
  itel <- 1
  repeat {
    anew <- template * crossprod(x, yold)
    svxa <- svd(x %*% anew)
    ynew <- tcrossprod(svxa$u[, 1:m], svxa$v[, 1:m])
    resi <- x - tcrossprod(ynew, anew)
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
}
