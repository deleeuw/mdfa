mdfaAlgorithmA <- function(x,
                           told,
                           proj = mdfaCFAProject,
                           itmax = 1000,
                           eps = 1e-10,
                           verbose = TRUE) {
  m <- nrow(told)
  svxa <- svd(x %*% told)
  yold <- tcrossprod(svxa$u[, 1:m], svxa$v[, 1:m])
  resi <- x - tcrossprod(yold, told)
  fold <- sum(resi^2)
  itel <- 1
  repeat {
    tnew <- proj(crossprod(x, yold))
    svxa <- svd(x %*% tnew)
    ynew <- tcrossprod(svxa$u[, 1:m], svxa$v[, 1:m])
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
  return(
    list(
      ydet = ynew,
      tmat = tnew,
      loss = fnew,
      itel = itel
    ))
}


