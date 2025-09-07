macdonell <- matrix(0, 7, 7)
macdonell[2, 1] <- .40163
macdonell[3, 1:2] <- c(.39454, .61779)
macdonell[4, 1:3] <- c(.30071, .15040, .32097)
macdonell[5, 1:4] <- c(.30539, .13515, .28869, .84638)
macdonell[6, 1:5] <- c(.33886, .20614, .36322, .75871, .79699)
macdonell[7, 1:6] <- c(.33993, .18308, .34527, .66084, .79986, .73636)
macdonell <- macdonell + t(macdonell)
diag(macdonell) <- 1
row.names(macdonell) <- c("Head Length", "Head Breadth", "Face Breadth", "Finger", "Cubit", "Foot", "Height")

mcaeig <- eigen(macdonell)
mcavec <- mcaeig$vectors[, 1:2]
mcaval <- diag(sqrt(mcaeig$values[1:2]))
mcacom <- mcavec %*% mcaval
mctold <- cbind(mcacom, diag(sqrt(1 - rowSums(mcacom^2))))

mcwdia <- diag(1 / diag(macdonell))
mcwful <- solve(macdonell)
mcwide <- diag(7)