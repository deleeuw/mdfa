library(numDeriv)

data(cattell, package = "psych")
dold <- 1 / diag(solve(cattell))

lsfa <- function(dold,
                 cmat,
                 p,
                 itmax = 5,
                 ieps = 1e-10,
                 iverbose = TRUE,
                 jtmax = 5,
                 jeps = 1e-15,
                 jverbose = TRUE) {
  itel <- 1
  eold <- eigen(cmat - diag(dold))
  fold <- sum(eold$values[-(1:p)]^2)
  aold <- eold$vectors[, 1:p] %*% sqrt(diag(eold$values[1:p]))
  repeat {
    dnew <- diag(cmat - tcrossprod(aold))
    enew <- eigen(cmat - diag(dnew))
    fnew <- sum(enew$values[-(1:p)]^2)
    anew <- enew$vectors[, 1:p] %*% sqrt(diag(enew$values[1:p]))
    if (iverbose) {
      cat(
        "itel ",
        formatC(itel, width = 4, format = "d"),
        "fold ",
        formatC(fold, digits = 15, format = "f"),
        "fnew ",
        formatC(fnew, digits = 15, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((fold - fnew) < ieps)) {
      break
    }
    itel <- itel + 1
    fold <- fnew
    aold <- anew
  }
  jtel <- 1
  dold <- dnew
  fold <- lsfaFunction(dold, cmat, p)
  repeat {
    g <- lsfaGradient(dold, cmat, p)
    h <- lsfaHessian(dold, cmat, p)
    dnew <- dold - solve(h, g)
    fnew <- lsfaFunction(dnew, cmat, p)
    if (jverbose) {
      cat(
        "jtel ",
        formatC(jtel, width = 4, format = "d"),
        "fold ",
        formatC(fold, digits = 15, format = "f"),
        "fnew ",
        formatC(fnew, digits = 15, format = "f"),"\n"
      )
    }
    if ((jtel == jtmax) || ((fold - fnew) < jeps)) {
      break
    }
    jtel <- jtel + 1
    dold <- dnew
    fold <- fnew
  }
  return(list(a = anew, d = dnew, f = fnew, itel = itel, jtel = jtel))
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
    for (eta in 1:n) {
      if (eta == nu) {
        next
      }
      lj <- lbd[eta]
      h <- h - 4 * aux * (1 / (lj - li)) * outer(vec[, eta], vec[, eta])
    }
  }
  return(h)
}

lsfaFunction <- function(theta, cmat, p) {
  return(sum(eigen(cmat - diag(theta))$values[-(1:p)]^2))
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
