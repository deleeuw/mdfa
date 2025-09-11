
data(bfi, package = "psych")
bfi <- as.matrix(bfi)[, 1:25]
n <- nrow(bfi)
m <- ncol(bfi)
q <- 5
p <- m + q
bfi <- e1071::impute(bfi, "mean")
bfi <- apply(bfi, 2, function(x) x - mean(x))
bfi <- apply(bfi, 2, function(x) x / sqrt(sum(x ^ 2)))
bfcor<- crossprod(bfi)
ceig <- eigen(cfi)
cvec <- ceig$vectors[, 1:q]
cval <- diag(sqrt(ceig$values[1:q]))
bfaold <- cvec %*% cval
bfdold <- diag(sqrt(1 - rowSums(bfaold^2)))
bftold <- cbind(bfaold, bfdold)
bfyold <- matrix(rnorm(n * p), n, p)
bfyold <- apply(bfyold, 2, function(x) x - mean(x))
bfyold <- qr.Q(qr(bfyold))
