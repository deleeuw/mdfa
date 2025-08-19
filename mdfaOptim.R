matrixPrint <- function(x,
                        digits = 6,
                        width = 8,
                        format = "f",
                        flag = "+") {
  print(noquote(
    formatC(
      x,
      digits = digits,
      width = width,
      format = format,
      flag = flag
    )
  ))
}

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
  
aeig <- eigen(emmett)
avec <- aeig$vectors[, 1:3]
aval <- diag(sqrt(aeig$values[1:3]))
acom <- avec %*% aval
aemm <- cbind(acom, diag(sqrt(1 - rowSums(acom^2))))
temm <- cbind(matrix(1, 9, 3), diag(9))

meig <- eigen(maxwell)
mvec <- meig$vectors[, 1:4]
mval <- diag(sqrt(meig$values[1:4]))
mcom <- mvec %*% mval
memm <- cbind(mcom, diag(sqrt(1 - rowSums(mcom^2))))
uemm <- cbind(matrix(1, 10, 4), diag(10))

optProj <- function(c, a, t) {
  nrow <- nrow(a)
  ncol <- ncol(a)
  ncom <- ncol - nrow
  h <- optim(
    par = as.vector(a),
    fn = funProj,
    gr = grdProj,
    method = "BFGS",
    control = list(trace = 6),
    cp = c,
    tp = t,
    nrow = nrow,
    ncol = ncol
  )
  hmat <- matrix(h$par, nrow, ncol)
  loadings <- hmat[, 1:ncom]
  uniquenesses <- diag(hmat[, -(1:ncom)]^2)
  mat <- crossprod(loadings, diag(1 / uniquenesses) %*% loadings)
  mvc <- eigen(mat)$vectors
  loadings <- loadings %*% mvc
  return(
    list(
      loadings = loadings,
      uniquenesses = uniquenesses,
      loss = h$value,
      counts = h$counts
    )
  )
}

funProj <- function(a, cp, tp, nrow, ncol) {
  tcp <- sum(diag(cp))
  tap <- sum(a^2)
  amt <- matrix(a, nrow, ncol)
  aca <- crossprod(amt, cp %*% amt)
  func <- tcp + tap - 2 * sum(sqrt(svd(aca)$d[1:nrow]))
  return(func)
}

grdProj <- function(a, cp, tp, nrow, ncol) {
  amt <- matrix(a, nrow, ncol)
  aca <- crossprod(amt, cp %*% amt)
  eca <- eigen(aca)
  evc <- eca$vectors[, 1:nrow]
  eva <- diag(1 / sqrt(eca$values[1:nrow]))
  mva <- evc %*% eva %*% t(evc)
  return(as.vector(tp * (a - cp %*% amt %*% mva)))
}

mdfaAlgorithmB <- function(c,
                    aold,
                    template,
                    itmax = 100,
                    eps = 1e-10,
                    verbose = TRUE) {
  itel <- 1
  m <- nrow(c)
  p <- ncol(aold)
  q <- p - m
  ssqc <- sum(diag(c))
  ssqa <- sum(aold^2)
  aca <- crossprod(aold, c %*% aold)
  eca <- eigen(aca)
  fold <- ssqc + ssqa - 2 * sum(sqrt(eca$values[1:m]))
  repeat {
    evc <- eca$vectors[, 1:m]
    eva <- sqrt(eca$values[1:m])
    mva <- evc %*% diag(1 / eva) %*% t(evc)
    anew <- template * (c %*% aold %*% mva)
    aca <- crossprod(anew, c %*% anew)
    eca <- eigen(aca)
    ssqa <- sum(anew^2)
    fnew <- ssqc + ssqa - 2 * sum(sqrt(eca$values[1:m]))
    if (verbose) {
      cat(
        "itel",
        formatC(itel, format = "d", width = 4),
        "fold",
        formatC(fold, digits = 10, format = "f"),
        "fnew",
        formatC(fnew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((fold - fnew) < eps)) {
      break
    }
    itel <- itel + 1
    aold <- anew
    fold <- fnew
  }
  loadings <- anew[, 1:q]
  uniquenesses <- diag(anew[, -(1:q)]) ^ 2
  mat <- crossprod(loadings, diag(1 / uniquenesses) %*% loadings)
  mvc <- eigen(mat)$vectors
  loadings <- loadings %*% mvc
  return(list(
    loadings = loadings,
    uniquenesses = uniquenesses,
    loss = fnew,
    itel = itel
  ))
}

mdfaAlgorithmA <- function(x,
                    aold,
                    template,
                    itmax = 1000,
                    eps = 1e-10,
                    verbose = TRUE) {
  m <- ncol(x)
  svxa <- svd(x %*% aold)
  yold <- tcrossprod(svxa$u[, 1:m], svxa$v[, 1:m])
  resi <- x - tcrossprod(yold, aold)
  fold <- sum(resi^2)
  itel <- 1
  repeat {
    anew <- template * crossprod(x, yold)
    svxa <- svd(x %*% anew)
    ynew <- tcrossprod(svxa$u[, 1:m], svxa$v[, 1:m])
    resi <- x - tcrossprod(ynew, anew)
    fnew <- sum(resi^2)
    if (verbose) {
      cat(
        "itel",
        formatC(itel, format = "d", width = 4),
        "fold",
        formatC(fold, digits = 10, format = "f"),
        "fnew",
        formatC(fnew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((fold - fnew) < eps)) {
      break
    }
    itel <- itel + 1
    yold <- ynew
    fold <- fnew
  }
}

matrixPower <- function(x, p = 1 / 2) {
  e <- eigen(x)
  evec <- e$vectors
  eval <- e$values
  epow <- eval^p
  return(evec %*% diag(epow) %*% t(evec))
}
