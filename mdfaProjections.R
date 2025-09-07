mdfaCFAProjection <- function(t) {
  q <- ncol(t) - nrow(t)
  return(cbind(t[, 1:q], diag(diag(t[, -(1:q)]))))
}

mdfaCFAWProjection <- function(t, wmat) {
  q <- ncol(t) - nrow(t)
  a <- t[, 1:q]
  dtilde <- diag(wmat %*% t[, -(1:q)]) / diag(wmat)
  return(cbind(a, diag(dtilde)))
}

cattellProjection <- function(tmat) {
  tmat[-(1:2), 2] <- 0
  tmat[-(3:4), 3] <- 0
  tmat[-(5:6), 4] <- 0
  tmat[-(7:8), 5] <- 0
  tmat[-(9:12), 6] <- 0
  tmat[, 7:18] <- diag(diag(tmat[, 7:18]))
  return(tmat)
}

cattellSVDProjection <- function(tmat) {
  tmat[1:2, 2:6] <- rankOne(tmat[1:2, 2:6])
  tmat[3:4, 2:6] <- rankOne(tmat[3:4, 2:6])
  tmat[5:6, 2:6] <- rankOne(tmat[5:6, 2:6])
  tmat[7:8, 2:6] <- rankOne(tmat[7:8, 2:6])
  tmat[9:12, 2:6] <- rankOne(tmat[9:12, 2:6])
  tmat[, 7:18] <- diag(diag(tmat[, 7:18]))
  return(tmat)
}

mdfaGroupCFAProjection <- function(tmat, ngfac, groups) {
  ngrps <- max(groups)
  ncomm <- ngfac + ngrps
  m <- nrow(tmat)
  for (k in 1:ngrps) {
    tmat[-which(groups == k), ngfac + k] <- 0
  }
  tmat[, (ncomm + 1):(ncomm + m)] <- diag(diag(tmat[, (ncomm + 1):(ncomm + m)]))
  return(tmat)
}
mdfaGroupCFACorFacProjection <- function () {}