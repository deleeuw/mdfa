

maxwell <- matrix(0, 10, 10)
maxwell[1, 2:10] <- c(.345, .594, .404, .579, .280, .449, .188, .303, .200)
maxwell[2, 3:10] <- c(.477, .338, .230, .159, .205, .120, .168, .145)
maxwell[3, 4:10] <- c(.498, .505, .251, .377, .186, .273, .154)
maxwell[4, 5:10] <- c(.389, .168, .249, .173, .195, .055)
maxwell[5, 6:10] <- c(.151, .285, .129, .159, .079)
maxwell[6, 7:10] <- c(.363, .359, .227, .260)
maxwell[7, 8:10] <- c(.448, .439, .511)
maxwell[8, 9:10] <- c(.429, .316)
maxwell[9, 10:10] <- c(.301)
maxwell <- maxwell + t(maxwell)
diag(maxwell) <- 1


meig <- eigen(maxwell)
mvec <- meig$vectors[, 1:4]
mval <- diag(sqrt(meig$values[1:4]))
mcom <- mvec %*% mval
memm <- cbind(mcom, diag(sqrt(1 - rowSums(mcom^2))))
uemm <- cbind(matrix(1, 10, 4), diag(10))
