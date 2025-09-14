set.seed(12345)
n <- 100
m <- 10
q <- 4
p <- 14
y <- matrix(rnorm(n * p), n, p)
y <- apply(y, 2, function(x) x - mean(x))
y <- qr.Q(qr(y))
f <- y[, 1:q]
u <- y[, -(1:q)]
a <- matrix(rnorm(m * q), m, q)
d <- diag(rnorm(m))
t <- cbind(a, d)
x <- tcrossprod(f, a) + u %*% d
arcov <- crossprod(x)
arcor <- cor(x)
ceig <- eigen(arcov)
cvec <- ceig$vectors[, 1:q]
cval <- diag(sqrt(ceig$values[1:q]))
araold <- cvec %*% cval
ardold <- diag(sqrt(diag(arcov) - rowSums(araold^2)))
artold <- cbind(araold, ardold)
aryold <- matrix(rnorm(n * p), n, p)
aryold <- apply(aryold, 2, function(x) x - mean(x))
aryold <- qr.Q(qr(aryold))
artemp <- cbind(matrix(1, m, q), diag(m))
