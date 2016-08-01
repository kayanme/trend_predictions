




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


form.draw <- function(x1, x2, x3, x4, y1, y2, y3, y4, var1, var2, var3, var4, x.borders, y.borders, shifted = F) {

    x <- seq(x.borders[1], x.borders[2], by = (x.borders[2] - x.borders[1]) / 1000)
    y <- seq(y.borders[1], y.borders[2], by = (y.borders[2] - y.borders[1]) / 1000)
    J <- y3 * x1 - y1 * x3 + y1 * x4 + y2 * x3 - y3 * x2 - y4 * x1 - y2 * x4 + y4 * x2
    print(J)  

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
        const.shift <- ((x3 - x4) * (y1 * x2 - y2 * x1) - (x1 - x2) * (y3 * x4 - y4 * x3))
        var.shift<-ya * (x1 - x2) - xa * (x3 - x4)             
        o.shift <- var.shift - const.shift
        o.shift <- o.shift / J

        o1 <- sqrt(var1 ^ 2 * (o.shift - x2) ^ 2 + var2 ^ 2 * (o.shift - x1) ^ 2)
        o2 <- sqrt(var3 ^ 2 * (o.shift - x4) ^ 2 + var4 ^ 2 * (o.shift - x3) ^ 2)
        print(c(const.shift = const.shift,
                mean.shift = mean(var.shift), sd.shift = sd(var.shift),
                mean.o1 = mean(o1), mean.o2 = mean(o2), sd.o1 = sd(o1), sd.o2 = sd(o2),
                range.o1=range(o1),range.o2=range(o2)))
        
        c1 <- exp( - (xa ^ 2 / (2 * o1 ^ 2) + ya ^ 2 / (2 * o2 ^ 2))) / (2 * pi * o1 * o2)

        c1
    }

    cut <- function(xa, ya) {
        o.shift <- - ((x3 - x4) * (y1 * x2 - y2 * x1) - (x1 - x2) * (y3 * x4 - y4 * x3))
        o.shift <- (ya * (x1 - x2) - xa * (x3 - x4) + o.shift) / J
        
        o1 <- sqrt(var1 ^ 2 * (o.shift - x2) ^ 2 + var2 ^ 2 * (o.shift - x1) ^ 2)
        o2 <- sqrt(var3 ^ 2 * (o.shift - x4) ^ 2 + var4 ^ 2 * (o.shift - x3) ^ 2)
       
        c1 <- (xa ^ 2 / (2*o1 ^ 2) + ya ^ 2 / (2*o2 ^ 2))
        c1*5
    }
    cut <- outer(x, y, FUN = cut)
    if (shifted) {
        e <- outer(x, y, FUN = f4)
    } else {
        e <- outer(x, y, FUN = f3)
    }
    plot3d(f4,xlim=c(-50,50),ylim=c(-50,50))
 #     contour(x, y, e, add = T, levels=c(0.0005,0.005,.01,.03,4), col = "red")
 #   contour(x, y, e, add = T, nlevels = 5, col = "red")
    #  contour(x, y, cut, add = T, levels = c(1, 20, 30,40, 50), col = "blue")
#    contour(x, y, cut, add = T, nlevels = 5, col = "blue")
}


pair.test <- function(a1, a2, a3, a4, x.borders = c(0, 1), y.borders = c(0, 1)) {
  
    par(mfrow = c(1, 1), mar = c(1, 1, 1, 1))
 
    plot.new()
    plot.window(xlim = x.borders, ylim = y.borders)
    test.data <- prepare.checks(a1[1:2], a1[3], a2[1:2], a2[3], a3[1:2], a3[3], a4[1:2], a4[3], test.count = 100000)
  
    cross.test.shifted.axis(a1[1:2], a1[3], a2[1:2], a2[3], a3[1:2], a3[3], a4[1:2], a4[3], test.data)
    cross.test.shifted(a1[1:2], a1[3], a2[1:2], a2[3], a3[1:2], a3[3], a4[1:2], a4[3], test.data)
    form.draw(a1[1], a2[1], a3[1], a4[1], a1[2], a2[2], a3[2], a4[2], a1[3], a2[3], a3[3], a4[3],
       x.borders = x.borders, y.borders = y.borders,shifted = T)

    plot.new()
    plot.window(xlim = x.borders, ylim = y.borders)
    cross.test.real.axis(a1[1:2], a1[3], a2[1:2], a2[3], a3[1:2], a3[3], a4[1:2], a4[3], test.data)
    cross.test.real(a1[1:2], a1[3], a2[1:2], a2[3], a3[1:2], a3[3], a4[1:2], a4[3], test.data)
    form.draw(a1[1], a2[1], a3[1], a4[1], a1[2], a2[2], a3[2], a4[2], a1[3], a2[3], a3[3], a3[3],
       x.borders = x.borders, y.borders = y.borders)
}

#pair.test(c(-2,-2,1),c(10,1,1),c(-5,5,1),c(0,0.05,.1),c(-100,100),c(-100,100))

pair.test(c(-2, -2, 1), c(10, 10, 10), c(-5, 5, 10), c(0, 0.05, .1), c(-100, 100), c(-100, 100))