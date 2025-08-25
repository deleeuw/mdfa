data(Tucker, package = "psych")

tuckerProjection <- function(tmat) {
  tmat[-(1:4), 3] <- 0
  tmat[-(5:9), 4] <- 0
  tmat[, 5:13] <- diag(diag(tmat[, 5:13]))
  return(tmat)
}

cmat <- as.matrix(Tucker)
aeig <- eigen(cmat)
avec <- aeig$vectors[, 1:4]
aval <- diag(sqrt(aeig$values[1:4]))
acom <- avec %*% aval
tutold <- cbind(acom, diag(sqrt(1 - rowSums(acom^2))))
tutold <- tuckerProjection(tutold)