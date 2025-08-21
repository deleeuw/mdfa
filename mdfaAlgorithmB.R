mdfaAlgorithmB <- function(c,
                           aold,
                           template,
                           itmax = 100,
                           eps = 1e-10,
                           verbose = TRUE) {
  itel <- 1
  m <- nrow(c)
  p <- ncol(aold)
  q <- p - m
  ssqc <- sum(diag(c))
  ssqa <- sum(aold^2)
  aca <- crossprod(aold, c %*% aold)
  eca <- eigen(aca)
  fold <- ssqc + ssqa - 2 * sum(sqrt(eca$values[1:m]))
  repeat {
    evc <- eca$vectors[, 1:m]
    eva <- sqrt(eca$values[1:m])
    mva <- evc %*% diag(1 / eva) %*% t(evc)
    anew <- template * (c %*% aold %*% mva)
    aca <- crossprod(anew, c %*% anew)
    eca <- eigen(aca)
    ssqa <- sum(anew^2)
    fnew <- ssqc + ssqa - 2 * sum(sqrt(eca$values[1:m]))
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
    aold <- anew
    fold <- fnew
  }
  loadings <- anew[, 1:q]
  uniquenesses <- diag(anew[, -(1:q)]) ^ 2
  mat <- crossprod(loadings, diag(1 / uniquenesses) %*% loadings)
  mvc <- eigen(mat)$vectors
  loadings <- loadings %*% mvc
  return(list(
    loadings = loadings,
    uniquenesses = uniquenesses,
    loss = fnew,
    itel = itel
  ))
}

