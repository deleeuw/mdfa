matrixPrint <- function(x,
                        digits = 6,
                        width = 8,
                        format = "f",
                        flag = "+") {
  print(noquote(
    formatC(
      x,
      digits = digits,
      width = width,
      format = format,
      flag = flag
    )
  ))
}

matrixPower <- function(x, p = 1 / 2) {
  e <- eigen(x)
  evec <- e$vectors
  eval <- e$values
  epow <- eval^p
  return(evec %*% diag(epow) %*% t(evec))
}

leftNullSpace <- function(x) {
  q <- qr(x)
  indi <- if (q$rank == 0L) {
    seq_len(ncol(x)) }
  else {
    -seq_len(q$rank)
  }
  return(qr.Q(q, complete = TRUE)[, indi, drop=FALSE])
}

mdfaConvertTtoAD <- function(tmat, rotate = TRUE) {
  q <- ncol(tmat) - nrow(tmat)
  loadings <- tmat[, 1:q]
  uniquenesses <- diag(tmat[, -(1:q)])^2
  if (rotate) {
    mat <- crossprod(loadings, diag(1 / uniquenesses) %*% loadings)
    mvc <- eigen(mat)$vectors
    loadings <- loadings %*% mvc
  }
  return(list(loadings = loadings, uniquenesses = uniquenesses))
}

mdfaCompleteY <- function(xmat, ymat, tmat) {
  m <- ncol(xmat)
  p <- ncol(ymat)
  q <- p - m
  hmat <- xmat %*% tmat
  kperp <- leftNullSpace(hmat)
  lperp <- leftNullSpace(t(hmat))
  return(ymat + tcrossprod(kperp[, 1:q], lperp))
}


