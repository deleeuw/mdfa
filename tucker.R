data(Tucker, package = "psych")

tuckerProjection <- function(tmat) {
  tmat[-(1:4), 3] <- 0
  tmat[-(5:9), 4] <- 0
  tmat[, 5:13] <- diag(diag(tmat[, 5:13]))
  return(tmat)
}

tugroups <- c(1, 1, 1, 1, 2, 2, 2, 2)

tucker <- as.matrix(Tucker)
aeig <- eigen(tucker)
avec <- aeig$vectors[, 1:4]
aval <- diag(sqrt(aeig$values[1:4]))
acom <- avec %*% aval
tutold <- cbind(acom, diag(sqrt(1 - rowSums(acom^2))))
tutold <- tuckerProjection(tutold)