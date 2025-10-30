library(numDeriv)
library(RSpectra)

source("data/emmett.R")
cmat <- emmett
m <- nrow(cmat)
p <- 3
dold <- (1 - p / (2 * m)) * (1 / diag(solve(cmat)))
hmat <- matrixPower(cmat, -1 / 2)


matrixPower <- function(x, p) {
  h <- eigen(x)
  hvec <- h$vectors
  hval <- diag(h$values^p)
  return(tcrossprod(hvec %*% hval, hvec))
}

f1 <- function(x) {
  return((1 / x) + log(x) - 1)
}

df1 <- function(x) {
  return(-1 / (x^2) + (1 / x))
}

ddf1 <- function(x) {
  return(2 / (x^3) - (1 / (x^2)))
}

ff1 <- list(f1, df1, ddf1)

f2 <- function(x) {
  return(0.5 * (1 - x)^2)
}

df2 <- function(x) {
  return(-(1 - x))
}

ddf2 <- function(x) {
  return(1)
}

ff2 <- list(f2, df2, ddf2)

f3 <- function(x) {
  return(0.5 * (log(x))^2)
}

df3 <- function(x) {
  return(log(x) / x)
}

ddf3 <- function(x) {
  return((1 - log(x)) / (x^2))
}

ff3 <- list(f3, df3, ddf3)

f4 <- function(x) {
  return(0.5 * ((x - 1)^2) / (x^2))
}

df4 <- function(x) {
  return((x - 1) / (x^3))
}

ddf4 <- function(x) {
  return((3 - (2 * x)) / x^4)
}

ff4 <- list(f4, df4, ddf4)

swainDerivatives <- function(theta, cmat, ff, p) {
  n <- length(theta)
  cmat <- matrixPower(cmat, -1 / 2)
  h <- eigen(cmat %*% diag(theta) %*% cmat)
  hval <- h$values
  hvec <- h$vectors
  y <- cmat %*% hvec
  g <- drop((y[, -(1:p)]^2) %*% ff[[2]](hval[-(1:p)]))
  h <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      for (nu in (p + 1):n) {
        h[i, j] <- h[i, j] + ff[[3]](hval[nu]) * (y[i, nu]^2) * (y[j, nu]^2)
        for (eta in 1:n) {
          if (eta == nu) {
            next
          }
          h[i, j] <- h[i, j] - 2 * ff[[2]](hval[nu]) / (hval[eta] - hval[nu]) * y[i, nu] * y[i, eta] * y[j, nu] * y[j, eta]
        }
      }
    }
  }
  return(list(g = g, h = h))
}

swainFunction <- function(theta, cmat, ff, p) {
  hmat <- matrixPower(cmat, -1 / 2)
  h <- eigen(hmat %*% diag(theta) %*% hmat)
  f <- sum(ff[[1]](h$values[-(1:p)]))
  return(f)
}

swainGradient <- function(theta, cmat, ff, p) {
  cmat <- matrixPower(cmat, -1 / 2)
  h <- eigen(cmat %*% diag(theta) %*% cmat)
  hval <- h$values
  hvec <- h$vectors
  y <- cmat %*% hvec
  g <- drop((y[, -(1:p)]^2) %*% ff[[2]](hval[-(1:p)]))
  return(g)
}

swainHessian <- function(theta, cmat, ff, p) {
  n <- length(theta)
  cmat <- matrixPower(cmat, -1 / 2)
  h <- eigen(cmat %*% diag(theta) %*% cmat)
  hval <- h$values
  hvec <- h$vectors
  y <- cmat %*% hvec
  h <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      for (nu in (p + 1):n) {
        h[i, j] <- h[i, j] + ff[[3]](hval[nu]) * (y[i, nu]^2) * (y[j, nu]^2)
        for (eta in 1:n) {
          if (eta == nu) {
            next
          }
          h[i, j] <- h[i, j] - 2 * ff[[2]](hval[nu]) / (hval[eta] - hval[nu]) * y[i, nu] * y[i, eta] * y[j, nu] * y[j, eta]
        }
      }
    }
  }
  return(h)
}

swainNumDerivatives <- function(theta, cmat, ff, p) {
  hmat <- matrixPower(cmat, -1 / 2)
  theFunc <- function(theta) {
    return(sum(ff[[1]](eigen(
      hmat %*% diag(theta) %*% hmat
    )$values)[-(1:p)]))
  }
  g <- jacobian(theFunc, theta)
  h <- hessian(theFunc, theta)
  return(list(g = g, h = h))
}

swain <- function(dold,
                  cmat,
                  ff,
                  p,
                  itmax = 5,
                  ieps = 1e-10,
                  iverbose = TRUE,
                  jtmax = 5,
                  jeps = 1e-15,
                  jverbose = TRUE) {
  itel <- 1
  eold <- eigs_sym(cmat - diag(dold), p)
  eval <- eold$values
  fold <- sum((cmat - diag(dold))^2) - sum(eval^2)
  aold <- eold$vectors %*% sqrt(diag(eval))
  repeat {
    dnew <- diag(cmat - tcrossprod(aold))
    enew <- eigs_sym(cmat - diag(dnew), p)
    eval <- enew$values
    fnew <- sum((cmat - diag(dnew))^2) - sum(eval^2)
    anew <- enew$vectors %*% sqrt(diag(eval))
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
  fold <- swainFunction(dold, cmat, ff, p)
  repeat {
    g <- swainGradient(dold, cmat, ff, p)
    h <- swainHessian(dold, cmat, ff, p)
    dnew <- dold - solve(h, g)
    fnew <- swainFunction(dnew, cmat, ff, p)
    if (jverbose) {
      cat(
        "jtel ",
        formatC(jtel, width = 4, format = "d"),
        "fold ",
        formatC(fold, digits = 15, format = "f"),
        "fnew ",
        formatC(fnew, digits = 15, format = "f"),
        "\n"
      )
    }
    if ((jtel == jtmax) || ((fold - fnew) < jeps)) {
      break
    }
    jtel <- jtel + 1
    dold <- dnew
    fold <- fnew
  }
  enew <- eigen(cmat - diag(dnew))
  anew <- enew$vectors[, 1:p] %*% diag(sqrt(enew$values[1:p]))
  mat <- crossprod(anew, diag(1 / dnew) %*% anew)
  mvc <- eigen(mat)$vectors
  anew <- anew %*% mvc
  return(list(
    a = anew,
    d = dnew,
    f = fnew,
    itel = itel,
    jtel = jtel
  ))
}

mySwainDerivatives <- function(theta, hmat, ff, p) {
  m <- nrow(hmat)
  eign <- eigen(hmat %*% diag(theta) %*% hmat)
  eval <- eign$values
  umat <- hmat %*% eign$vectors
  alph <- ff[[2]](eval)
  beta <- 2 * alph * (1 / (outer(eval, eval, "-") + diag(m)))
  diag(beta) <- ff[[3]](eval)
  grad <- drop((umat[, -(1:p)]^2) %*% alph[-(1:p)])
  hess <- matrix(0, m, m)
  for (r in 1:m) {
    for (s in 1:m) {
      for (i in (p + 1):m) {
        for (j in 1:m) {
          hess[r, s] <- hess[r, s] + beta[i, j] * umat[r, i] * umat[r, j] * umat[s, i] * umat[s, j]
        }
      }
    }
  }
  return(list(grad = grad, hess = hess))
}
