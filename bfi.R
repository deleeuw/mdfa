
data(bfi, package = "psych")
bfi <- as.matrix(bfi)[, 1:25]
bfi <- e1071::impute(bfi, "mean")
bfi <- apply(bfi, 2, function(x) x - mean(x))
bfi <- apply(bfi, 2, function(x) x / sqrt(sum(x ^ 2)))
cfi <- cor(bfi)
ceig <- eigen(cfi)
cvec <- ceig$vectors[, 1:5]
cval <- diag(sqrt(ceig$values[1:5]))
aold <- cvec %*% cval
dold <- diag(sqrt(1 - rowSums(ccom^2)))
told <- cbind(ccom, demm)
