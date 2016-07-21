
t <- function(x1, y1, x2, y2) {
    cbind(
      k = (y1 - y2) / (x1 - x2),
      b = -(x1 / (x1 - x2)) * (y1 - y2) + y1
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


form.draw <- function(x1, x2, x3, x4, y1, y2, y3, y4, var1, var2, var3, var4, x.borders, y.borders,shifted=F) {

    x <- seq(x.borders[1], x.borders[2], by = .01)
    y <- seq(y.borders[1], y.borders[2], by = .001)
    J <- y3 * x1 - y1 * x3 + y1 * x4 + y2 * x3 - y3 * x2 - y4 * x1 - y2 * x4 + y4 * x2
    print(J)
    #v1 <- c(
           #(x1 - x2),
           #(y1 - y2))
    #v2 <- c(
          #(x3 - x4),
          #(y3 - y4))
    #cross <- c(9.75, -4.625)
    #print(v1)
    #print(v2)
    #axis <- matrix(
        #c(cross, cross + v1 / J * 4,
        #cross, cross - v1 / J * 4,
        #cross, cross + v2 / J * 4,
        #cross, cross - v2 / J * 4),
        #ncol = 4, byrow = T)
    #segments(x0 = axis[, 1], y0 = axis[, 2], x1 = axis[, 3], y1 = axis[, 4], col = "blue")
    #print("pred axis")

    f3 <- function(xa, ya) {

        o1 <- sqrt(var1 ^ 2 * (xa - x2) ^ 2 + var2 ^ 2 * (xa - x1) ^ 2)
        o2 <- sqrt(var3 ^ 2 * (xa - x4) ^ 2 + var4 ^ 2 * (xa - x3) ^ 2)

        c1 <- exp( - (xa * (y1 - y2) + x1 * (y2 - ya) + x2 * (ya - y1)) ^ 2 / (2 * o1 ^ 2)
        - (xa * (y3 - y4) + x3 * (y4 - ya) + x4 * (ya - y3)) ^ 2 / (2 * o2 ^ 2)) / (2 * pi * o1 * o2)
        a <- cbind(xa, ya, o1, o2, c1)
        b <- a[c1 == max(c1), c(1, 2, 3, 4, 5)]

        c1
    }

    f4 <- function(xa, ya) {
        o.shift <- (x3 - x4) * (y1 * x2 - y2 * x1) - (x1 - x2) * (y3 * x4 - y4 * x3) - xa * (x3 - x4) - ya * (x1 - x2)

        o1 <- sqrt(var1 ^ 2 * (x2 * J + o.shift) ^ 2 + var2 ^ 2 * (x1 * J + o.shift) ^ 2) / J
        o2 <- sqrt(var3 ^ 2 * (x4 * J + o.shift) ^ 2 + var4 ^ 2 * (x3 * J + o.shift) ^ 2) / J

        c1 <- exp( - xa ^ 2 / (2 * o2 ^ 2)
                   - ya ^ 2 / (2 * o1 ^ 2)) / (2 * pi * o1 * o2)

        c1
    }
    if (shifted) {
        e <- outer(x, y, FUN = f4)
    }
        else {
            e <- outer(x, y, FUN = f3)
        }
    contour(x, y, e, add = T, levels=c(0.05,.1,.3,3,7.9), col = "red")

}


pair.test <- function(a1, a2, a3, a4, x.borders = c(0, 1), y.borders = c(0, 1)) {
  
    par(mfrow = c(1, 1), mar = c(1, 1, 1, 1))
 
    plot.new()
    plot.window(xlim = x.borders, ylim = y.borders)
    test.data <- prepare.checks(a1[1:2], a1[3], a2[1:2], a2[3], a3[1:2], a3[3], a4[1:2], a4[3], test.count = 100000)
#    cross.test.data(test.data)
    cross.test.shifted(a1[1:2], a1[3], a2[1:2], a2[3], a3[1:2], a3[3], a4[1:2], a4[3], test.data)
    form.draw(a1[1], a2[1], a3[1], a4[1], a1[2], a2[2], a3[2], a4[2], a1[3], a2[3], a3[3], a4[3],
       x.borders = x.borders, y.borders = y.borders,shifted = T)

    #plot.new()
    #plot.window(xlim = x.borders, ylim = y.borders)
    #cross.test.real(a1[1:2], a1[3], a2[1:2], a2[3], a3[1:2], a3[3], a4[1:2], a4[3], test.data)
    #form.draw(a1[1], a2[1], a3[1], a4[1], a1[2], a2[2], a3[2], a4[2], a1[3], a2[3], a3[3], a3[3],
       #x.borders = x.borders, y.borders = y.borders)
}

pair.test(c(-2,2,.01),c(-1,1,.01),c(-2,-2,.1),c(-1,-1,.1),c(-1,1),c(-.1,.1))