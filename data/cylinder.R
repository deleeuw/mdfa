cylinders <- matrix(
  c(
    1,
    2,
    1,
    2,
    2,
    1,
    3,
    2,
    1,
    1,
    3,
    1,
    2,
    3,
    1,
    3,
    3,
    1,
    1,
    4,
    1,
    2,
    4,
    1,
    3,
    4,
    1,
    1,
    2,
    2,
    2,
    2,
    2,
    3,
    2,
    2,
    1,
    3,
    2,
    2,
    3,
    2,
    3,
    3,
    2,
    1,
    4,
    2,
    2,
    4,
    2,
    3,
    4,
    2,
    1,
    2,
    3,
    2,
    2,
    3,
    3,
    2,
    3,
    1,
    3,
    3,
    2,
    3,
    3,
    3,
    3,
    3,
    1,
    4,
    3,
    2,
    4,
    3,
    3,
    4,
    3
  ),
  27,
  3,
  byrow = TRUE
)
colnames(cylinders) <- c("d", "l", "s")
cylinder <- cbind(
  cylinders[, 1],
  cylinders[, 2],
  pi * cylinders[, 1]^2 / 4,
  pi * cylinders[, 1] * cylinders[, 2],
  pi * cylinders[, 1]^2 * cylinders[, 2] / 4,
  sqrt(cylinders[, 1]^2 + cylinders[, 2]^2),
  cylinders[, 3] * pi * cylinders[, 1]^2 * cylinders[, 2] / 4
)
cylinder <- apply(cylinder, 2, function(x) x - mean(x))
cylinder <- apply(cylinder, 2, function(x) x / sqrt(sum(x^2)))
m <- 7
q <- 2
p <- 9
n <- 27
cycor <- crossprod(cylinder)
cyuni <- 1 - c(.978, 1.000, .987, .968, .962, .968, .672)
ceig <- eigen(cycor - diag(cyuni))
cvec <- ceig$vectors[, 1:q]
cval <- diag(sqrt(ceig$values[1:q]))
cyaold <- cvec %*% cval
cydold <- diag(sqrt(cyuni))
cydold <- .5 * diag(7)
cytold <- cbind(cyaold, cydold)


