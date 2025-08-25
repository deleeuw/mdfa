data(cattell, package = "psych")

cattellProjection <- function(tmat) {
  tmat[-(1:2), 2] <- 0
  tmat[-(3:4), 3] <- 0
  tmat[-(5:6), 4] <- 0
  tmat[-(7:8), 5] <- 0
  tmat[-(9:12), 6] <- 0
  tmat[, 7:18] <- diag(diag(tmat[, 7:18]))
  return(tmat)
}

aeig <- eigen(cattell)
avec <- aeig$vectors[, 1:6]
aval <- diag(sqrt(aeig$values[1:6]))
acom <- avec %*% aval
aemm <- cbind(acom, diag(sqrt(1 - rowSums(acom^2))))
ctmat <- cattellProjection(aemm)

