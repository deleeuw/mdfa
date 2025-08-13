library(numDeriv)

set.seed(12345)
n = 100
m = 4
p = 6
q = 2
a <- cbind(matrix(rnorm(m * q), m, q), diag(rnorm(m)))
avec <- as.vector(a)
x <- matrix(rnorm(m * n), n, m)
x <- apply(x, 2, function(x)
  x - mean(x))
c <- crossprod(x) / n
aca <- crossprod(a, c %*% a)
eca <- eigen(aca)
lbd <- sqrt(eca$values[1:m])
kve <- eca$vectors[, 1:m]
f <- sum(lbd)
aca2 <- kve %*% diag(1 / lbd) %*% t(kve)
df <- c %*% a %*% aca2
eps <- 1e-6
ef <- matrix(0, m * p, m * p)
k <- 1
for (t in 1:p) {
  for (l in 1:m) {
    delta <- matrix(0, m, p)
    delta[l, t] <- 1
    aux <- crossprod(delta, c %*% a)
    hnum <- crossprod(kve, (aux + t(aux)) %*% kve)
    hden <- outer(lbd, lbd) * outer(lbd, lbd, "+")
    h <- kve %*% (hnum / hden) %*% t(kve)
    ef[k, ] <- as.vector(c %*% delta %*% aca2 - c %*% a %*% h)
    k <- k + 1
  }
}

theFunc <- function(avec, cmat = c, arow = m , acol = p) {
  amat <- matrix(avec, arow, acol)
  aca <- crossprod(amat, cmat %*% amat)
  return(sum(sqrt(eigen(aca)$values[1:arow])))
}

numCheck <- function(a, c) {
  m <- nrow(a)
  p <- ncol(a)
  mp <- m * p
  avec <- as.vector(a)
  g <- grad(theFunc, avec, cmat = c, arow = m, acol = p)
  h <- hessian(theFunc, avec, cmat = c, arow = m, acol = p)
  return(list(grad = matrix(g, m, p), hessian = matrix(h, mp, mp)))
}
