library(numDeriv)
m <- 4
p <- 6
a <- cbind(matrix(rnorm(8), 4, 2), diag(rnorm(4)))
cc <- crossprod(matrix(rnorm(400), 100, 4))/100

iter <- function(a, cc, itmax = 100, eps = 1e-6, verbose = TRUE) {
  itel <- 1
  m <- nrow(a)
  p <- ncol(a)
  q <- p - m
  sc <- sum(cc^2)
  oloss <- Inf
  repeat {
    sa <- sum(a^2)
    aca <- t(a) %*% cc %*% a
    ece <- eigen(aca)
    eca <- ece$values[1:m]
    ecv <- ece$vectors[, 1:m]
    nloss <- sc + sa - 2 * sum(sqrt(eca))
    acb <- ecv %*% diag(1  / sqrt(eca)) %*% t(ecv)
    fd <- cc %*% a %*% acb
    gg <- max(abs(a - fd))
    if (verbose) {
      cat("itel ", formatC(itel, format = "d"),
          "oloss ", formatC(oloss, digits = 15, format = "f"),
          "nloss ", formatC(nloss, digits = 15, format = "f"),
          "grad ", formatC(gg, digits = 15, format = "f"), "\n")
    }
    if ((itel == itmax) || ((oloss - nloss) < eps)) {
      break
    }
    itel <- itel + 1
    oloss <- nloss
    a[, 1:q] <- fd[, 1:q]
    a[, q + 1:m] <- diag(diag(fd[, q + 1:m]))
  }
  return(list(a = a, loss = nloss, itel = itel))
}

numD <- function(a, cc) {
  ff <- function(x) {
    amat <- matrix(x, m, p)
    aca <- t(amat) %*% cc %*% amat
    ece <- eigen(aca)
    return(sum(sqrt(ece$values[1:m])))
  } 
  fn <- grad(func = ff, x = as.vector(a))
  return(matrix(fn, 4, 6))
}

perD <- function(a, cc, eps) {
  m <- nrow(a)
  p <- ncol(a)
  ca <- t(a) %*% cc %*% a
  ea <- eigen(ca)
  sa <- sum(sqrt(ea$values[1:4]))
  de <- matrix(0, m, p)
  for (j in 1:m) {
    for (s in 1:p) {
      ap <- a
      ap[j, s] <- ap[j, s] + eps
      cp <- t(ap) %*% cc %*% ap
      ep <- eigen(cp)
      sp <- sum(sqrt(ep$values[1:4]))
      de[j, s] <- (sp - sa) / eps
    }
  }
return(de)
}

