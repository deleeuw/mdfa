set.seed(12345)
n = 100
m = 4
p = 6
q = 2
a <- cbind(matrix(rnorm(m * q), m, q), diag(rnorm(m)))
x <- matrix(rnorm(m * n), n, m)
x <- apply(x, 2, function(x) x - mean(x))
c <- crossprod(x) / n
aca <- crossprod(a, c %*% a)
eca <- eigen(aca)
lbd <- sqrt(eca$values[1:m])
kve <- eca$vectors[, 1:m]
f <- sum(lbd)
aca2 <- kve %*% diag(1 / lbd) %*% t(kve)
df <- c %*% a %*% aca2
eps <- 1e-6
ef <- matrix(0, m, p)
for (j in 1:m) {
  for (s in 1:p) {
    del <- a
    del[j, s] <- del[j, s] + eps
    acap <- crossprod(del, c %*% del)
    ecap <- eigen(acap)
    lbdp <- sqrt(ecap$values[1:m])
    kvep <- ecap$vectors[, 1:m]
    fp <- sum(lbdp)
    ef[j, s] <- (fp - f) / eps
  }
}
