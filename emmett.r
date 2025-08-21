emmett <- matrix(0, 9, 9)
emmett[1, 2:9] <- c(0.523, 0.395, 0.471, 0.346, 0.426, 0.576, 0.434, 0.639)
emmett[2, 3:9] <- c(0.479, 0.506, 0.418, 0.462, 0.547, 0.283, 0.645)
emmett[3, 4:9] <- c(0.355, 0.270, 0.254, 0.452, 0.219, 0.504)
emmett[4, 5:9] <- c(0.691, 0.791, 0.443, 0.285, 0.505)
emmett[5, 6:9] <- c(0.679, 0.383, 0.149, 0.409)
emmett[6, 7:9] <- c(0.372, 0.314, 0.472)
emmett[7, 8:9] <- c(0.385, 0.680)
emmett[8, 9:9] <- c(0.470)
emmett <- emmett + t(emmett)
diag(emmett) <- 1

aeig <- eigen(emmett)
avec <- aeig$vectors[, 1:3]
aval <- diag(sqrt(aeig$values[1:3]))
acom <- avec %*% aval
aemm <- cbind(acom, diag(sqrt(1 - rowSums(acom^2))))
temm <- cbind(matrix(1, 9, 3), diag(9))
