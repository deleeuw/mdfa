library(fungible)
data(AmzBoxes)

rankMe <- function(x) {
  n <- length(x)
  return(ifelse(outer(1:n, order(x), "=="), 1, 0) %*% (1:n))
}

BoxList <- GenerateBoxData (XYZ = AmzBoxes[,2:4],
                            BoxStudy = 20,  
                            Reliability = .75,
                            ModApproxErrVar = .10,
                            SampleSize = NULL, 
                            NMinorFac = 50,
                            epsTKL = .20,
                            Seed = 1,
                            SeedErrorFactors = 1,
                            SeedMinorFactors = 2,
                            PRINT = FALSE,
                            LB = FALSE,
                            LBVal = 1,
                            Constant = 0)
box0 <- BoxList$BoxData[, 1:17]
box1 <- apply(box0, 2, function(x) x - mean(x))
box2 <- apply(box1, 2, function(x) {x / sqrt(sum(x^2))})
boxcv <- cov(box0)
boxco <- cor(box0)
ceig <- eigen(boxco)
q <- 3
cvec <- ceig$vectors[, 1:q]
cval <- diag(sqrt(ceig$values[1:q]))
aold <- cvec %*% cval
dold <- diag(sqrt(1 - rowSums(aold^2)))
bold <- cbind(aold, dold)
