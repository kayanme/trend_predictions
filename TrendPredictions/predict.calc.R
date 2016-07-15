
t <- function(x1, y1, x2, y2) {
    cbind(
      k = (y1 - y2) / (x1 - x2),
      b = (x1 / (x1 - x2)) * (y2 + y1) + y1
      )
}


t.var <- function(x1, var1, x2, var2) {
    cbind(
      k = (var1 + var2) / ((x1 - x2)^2),
      b = (x1 / (x1 - x2) +1) ^ 2 * var1 + (x1 / (x1 - x2)) ^ 2 * var2
      )

}

cross.cov <- function(varb1, vark1, varb2, vark2, r) {

}

t.cov <- function(x1,  x2,var1,var2) {
    e <- 1 / (x1 - x2)
    e*(e*x1+1)*var1 - e*e*x2*var2
}

test.sd <- function(c1, c2, sd1, sd2,x1,x2,test.count=100000) {
    y1 <- qnorm(runif(1:test.count, 0, 1), mean = c1, sd = sd1)
    y2 <- qnorm(runif(1:test.count, 0, 1), mean = c2, sd = sd2)
    coeffs <- t(x1, y1, x2,y2)
    fact.mean.k <- mean(coeffs[, "k"])
    fact.mean.b <- mean(coeffs[, "b"])
    fact.sd.k <- var(coeffs[, "k"])
    fact.sd.b <- var(coeffs[, "b"])

    pred <- t(x1, c1, x2, c2)
    pred2 <- t.var(x1, sd1^2, x2, sd2^2)

    pred.mean.k <- pred[1,"k"]
    pred.mean.b <- pred[1,"b"]
    pred.sd.k <- pred2[1,"k"]
    pred.sd.b <- pred2[1, "b"]
    c <- t.cov(x1, x2, sd1 ^ 2, sd2 ^ 2)# /  (sd1 * sd2)

    e<-matrix(c(fact.mean.k, fact.mean.b, fact.sd.k, fact.sd.b, pred.mean.k, pred.mean.b, pred.sd.k, pred.sd.b), ncol = 2)
    colnames(e) <- c("fact", "pred")
    rownames(e) <- c("mean.k", "mean.b", "var.k", "var.b")
    
    list(base = e, fact.cov = cor(coeffs), pred.cov = matrix(c(1,c,c,1),byrow=T,ncol=2))
}



cross <- function(k1, b1, k2, b2) {
    cbind(
          x = (b2 - b1) / (k1 - k2),
          y = (b2 - b1) / (k1 - k2)*k1+b1
           )

}

test.cross <- function(mean1, mean2, cov1, cov2,  test.count = 1000000) {
    f1.k <- qnorm(runif(1:test.count, 0, 1), mean = mean1["k"], sd = sqrt(cov1[1,1]))
    f1.b <- qnorm(runif(1:test.count, 0, 1), mean = mean1["b"], sd = sqrt(cov1[2, 2])) + f1.k * cov1[1, 2] 
    f1 <- cbind(f1.k, f1.b)

    f2.k <- qnorm(runif(1:test.count, 0, 1), mean = mean2["k"], sd = sqrt(cov2[1, 1]))
    f2.b <- qnorm(runif(1:test.count, 0, 1), mean = mean2["b"], sd = sqrt(cov2[2, 2])) + f2.k * cov2[1, 2]
    f2 <- cbind(f2.k, f2.b)

    coeffs <- cross(f1.k,f1.b,f2.k,f2.b)

    fact.mean.x <- mean(coeffs[, "x"])
    fact.mean.y <- mean(coeffs[, "y"])
    fact.cov <- cor(coeffs)

    pred <- cross(mean(f1.k), mean(f1.b), mean(f2.k), mean(f2.b))

    pred.mean.x <- pred[1, "x"]
    pred.mean.y <- pred[1, "y"]

    e <- matrix(c(fact.mean.x, fact.mean.y, pred.mean.x, pred.mean.y), ncol = 2)
    colnames(e) <- c("fact", "pred")
    rownames(e) <- c("x", "y")
    list(means=e,fact.cov=fact.cov)
}

print(test.sd(20,20,.01,.01,4,3,1000000))
#print(var(1/qnorm(runif(1:100000, 0, 1), mean = 1, sd = 3)))
#print(test.cross(c(k = 1, b = 0), c(k = -1, b = 5), rbind(c(.1, -.3), c(-.3, 1)), rbind(c(.1, -.3), c(-.3, 1))))

#k <- 2
#b <- 0
#r<-.99
#v.k <- 1
#v.b<- 50
#par(mfrow=c(1,1))
#x <- c(-1000:1000)
#mid.y <- k * x + b
#var.y <- sqrt((v.b ^ 2) - 2 * r * v.k * v.b * x + v.k ^ 2 * x ^ 2)


#plot(x, mid.y, type = "l",ylim=c(-100,100),xlim=c(-1000,1000))
#lines(x, mid.y + var.y, type = "l",col="red")
#lines(x, mid.y - var.y, type = "l", col = "green")


