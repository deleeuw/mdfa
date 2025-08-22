mdfaAlgorithmA <- function(x,
                           aold,
                           dold,
                           projA = mdfaProjectA,
                           projD = mdfaProjectD,
                           itmax = 1000,
                           eps = 1e-10,
                           verbose = TRUE) {
  m <- ncol(x)
  q <- ncol(aold)
  told <- cbind(aold, dold)
  svxa <- svd(x %*% told)
  yold <- tcrossprod(svxa$u[, 1:m], svxa$v[, 1:m])
  resi <- x - tcrossprod(yold, told)
  fold <- sum(resi^2)
  itel <- 1
  repeat {
    tnew <- crossprod(x, yold)
    anew <- projA(tnew[, 1:q])
    dnew <- projD(tnew[, -(1:q)])
    tnew <- cbind(anew, dnew)
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
      f = ynew[, 1:q],
      u = ynew[, -(1:q)],
      loadings = anew,
      uniquenesses = diag(dnew),
      loss = fnew,
      itel = itel
    ))
}

mdfaProjectA <- function(a) {
  return(a)
}

mdfaProjectD <- function(d) {
  return(diag(diag(d)))
}