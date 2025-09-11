set.seed(12345)
y <- matrix(rnorm(1400), 100, 14)
y <- apply(y, 2, function(x) x - mean(x))
y <- qr.Q(qr(y))
f <- y[, 1:4]
u <- y[, -(1:4)]
a <- matrix(rnorm(40), 10, 4)
d <- diag(rnorm(10))
t <- cbind(a, d)
x <- tcrossprod(f, a) + u %*% d
n <- 100
m <- 10
q <- 4
p < 14
arcov <- cov(x)
arcor <- cor(x)
ceig <- eigen(arcor)
cvec <- ceig$vectors[, 1:q]
cval <- diag(sqrt(ceig$values[1:q]))
araold <- cvec %*% cval
ardold <- diag(sqrt(1 - rowSums(araold^2)))
artold <- cbind(araold, ardold)
aryold <- matrix(rnorm(n * p), n, p)
aryold <- apply(aryold, 2, function(x) x - mean(x))
aryold <- qr.Q(qr(aryold))
