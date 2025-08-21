library(numDeriv)

data(cattell, package = "psych")

lsfaDerivatives <- function(theta, cmat, p) {
  n <- nrow(cmat)
  a <- cmat - diag(theta)
  h <- eigen(a)
  indi <- (p + 1):n
  lbd <- h$values
  vec <- h$vectors
  f <- sum(lbd[indi]^2)
  g <- -2 * drop(vec[, indi]^2 %*% lbd[indi])
  h <- 2 * tcrossprod(vec[, indi]^2)
  for (nu in indi) {
    li <- lbd[nu]
    aux <- li * outer(vec[, nu], vec[, nu])
    for(eta in 1:n) {
      if (eta == nu) {
        next
      }
      lj <- lbd[eta]
      h <- h - 4 * aux * (1 / (lj - li)) * outer(vec[, eta], vec[, eta])
    }
  }
  return(list(f = f, g = g, h = h))
}

lsfaGradient <- function(theta, cmat, p) {
  n <- length(theta)
  h <- eigen(cmat - diag(theta))
  indi <- (p + 1):n
  lbd <- h$values
  vec <- h$vectors
  g <- -2 * drop(vec[, indi]^2 %*% lbd[indi])
  return(g)
}

lsfaHessian <- function(theta, cmat, p) {
  n <- nrow(cmat)
  a <- cmat - diag(theta)
  h <- eigen(a)
  indi <- (p + 1):n
  lbd <- h$values
  vec <- h$vectors
  h <- 2 * tcrossprod(vec[, indi]^2)
  for (nu in indi) {
    li <- lbd[nu]
    aux <- li * outer(vec[, nu], vec[, nu])
    for(eta in 1:n) {
      if (eta == nu) {
        next
      }
      lj <- lbd[eta]
      h <- h - 4 * aux * (1 / (lj - li)) * outer(vec[, eta], vec[, eta])
    }
  }
  return(h)
}

lsfaNumDerivatives <- function(theta, cmat, p) {
  theFunc <- function(theta) {
    n <- length(theta)
    h <- eigen(cmat - diag(theta))
    indi <- (p + 1):n
    f <- sum(h$values[indi]^2)
    return(f)
  }
  f <- theFunc(theta)
  g <- grad(theFunc, theta)
  h <- hessian(theFunc, theta)
  return(list(f = f, g = g, h = h))
}
