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

matrixPower <- function(x, p = 1 / 2) {
  e <- eigen(x)
  evec <- e$vectors
  eval <- e$values
  epow <- eval^p
  return(evec %*% diag(epow) %*% t(evec))
}