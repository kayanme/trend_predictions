

x.range = c(1,1000)

noise.error<-.05
noise.amp<-1000

#data.trends<-c(
                    #1, -4,
                    #1000, 5,
                    #2000, -2,
                    #3000, 3,
                    #4000, -5,
                    #5000, -1,
                    #6000, 0,
                    #7000, 10,
                    #8000, 0,
                    #9000, -8)
#data.trends <- c(1,5)

data.trends <- c(
                    1, -4,                    
                    500, 5
                    )




x.values <- c(x.range[1]:x.range[2])
trends <- matrix(
                 data = data.trends,
                 by = T,ncol=2
)

colnames(trends)<-c("x","k")



#noises <- matrix(data=c(noise1,noise2,noise3,noise4),byrow = T,ncol=2)

trend.ext <- cbind(x = trends[, "x"], 
                   xl = c(trends[, "x"],max(x.values))[2:(nrow(trends)+1)],
                   k = trends[, "k"])
init.shift<-x.values[1]*trends[1,"k"]
trend.ext <- cbind(trend.ext, b = (trend.ext[, "xl"] - trend.ext[, "x"]) * trend.ext[, "k"] )
trend.ext <- cbind(x = trend.ext[, "x"],
                   xl = trend.ext[, "xl"],
                   k = trend.ext[, "k"],
                   b = apply(trend.ext, MARGIN = 1, FUN =
                                                 function(x) {
                                                     sum(trend.ext[trend.ext[, "x"] < x["x"], "b"]) - init.shift - x["x"] * x["k"]
                                                 }))


trend.lines <- sapply(x.values, FUN = function(x) {
    r <- trend.ext[trend.ext[, "x"] <= x & trend.ext[, "xl"] >= x,]
    if (length(r) > 4) {r<-r[1,] }
    x*r["k"] + r["b"]
})


noised <- qnorm(runif(length(x.values), 0,1))*trend.lines*noise.error +  trend.lines



ma <- function(arr, n) {
    halfRange <- (n - 1) / 2
    res = arr
    disp = arr
    arr <- c(seq(0, 0, length = halfRange), arr, seq(0, 0, length = halfRange))
 
 
    conv <- seq( - (n - 1) / 2, (n - 1) / 2) * (pi / (2 * n)) 
    conv <- cos(conv)*1.11
    conv<-seq(1,1,length(n))
    for (i in 1:length(res)) {
        part <- arr[i:(i+halfRange*2)]
        res[i] = mean(part * conv)
        disp[i] = sqrt(sum((part-res[i])^2)/(n-1))
    }
    
   cbind(mean = res,sd=disp) 
}



approx <- function(x, y) {
    xshift <- x[1]
    yshift<-y[1]
    x <- x - xshift
    y <- y - yshift
    x1 <- length(x)
    x2 <- sum(x)
    x3 <- sum(as.numeric(x) * as.numeric(x))
    y1 <- sum(y)
    y2 <- sum(y * x)

  
    t<-solve(
          matrix(c(x1, x2, x2, x3), ncol = 2, byrow = T),
          c(y1,y2)
          )   
    k <- t[2]
    b <- t[1] - k * xshift + yshift
    list(y=x*k+b,k=k,b=b)
}



where.point <- function(p1, p2, p3) {
    v <- c(p2[1] - p1[1], p2[2] - p1[2])

    out <- p2 + v * (p3[1] - p2[1])

    p3[2] - out[2]

}

are.crossed.right <- function(a1, a2, b1, b2) {
    (a1[1] - b1[1]) ^ 2 + (a1[2] - b1[2]) ^ 2 > (a2[1] - b2[1]) ^ 2 + (a2[2] - b2[2]) ^ 2
}

swinging.door <- function(x, n) {
    halfPart <- (n - 1) / 2
    means <- x[, "mean"]
    sds <- x[, "sd"]

    std.err <- sds[1] / sqrt(n)
    u1 <- c(1, means[1] + std.err)
    l1 <- c(1, means[1] - std.err)
    u2 <- c(2, means[2])
    l2 <- c(2, means[2])

    res <- c(1, x[1], std.err, 0)
    for (i in 1:(length(means))) {
        if (where.point(u1, u2, c(i, means[i])) > 0) {
            u2 <- c(i, means[i])
        }

        if (where.point(l1, l2, c(i, means[i])) < 0) {
            l2 <- c(i, means[i])
        }

        if (!are.crossed.right(u1, u2, l1, l2)) {

            trend <- res[length(res) - 3]:i
            trend.length <- length(trend)

            dev.seq <- if (trend.length > 2 * n) {
                trend
            } else {
                a <- (i - 2 * n)
                if (a > 0) {
                    a:i
                }
                    else {
                        1:i
                    }
            }

            std.err <- sqrt(mean(sds[dev.seq] ^ 2)) / sqrt(n)            
            final.value.appox <- approx(trend, means[trend])
            final.value <- final.value.appox$k * (i - 1) + final.value.appox$b

            if (abs(final.value.appox$k - res[length(res)]) > 0.1 & trend.length > halfPart * 2) {
                res <- c(res, i - 1, final.value, std.err, final.value.appox$k)
            }
            #if (trend.length > halfPart * 2) {
                #res <- c(res, i - 1, means[i - 1], std.err, final.value.appox$k)
            #}

            u1 <- c(i - 1, final.value + std.err )
            l1 <- c(i - 1, final.value - std.err )
            u2 <- c(i, means[i])
            l2 <- c(i, means[i])
        }
    }

    trend <- length(means):res[length(res) - 3]
    trend.length <- length(trend)

    std.err <- sqrt(sum(sds[trend] ^ 2)) / trend.length

    approxed <- approx(trend, means[trend])

    if (length(means) - res[length(res) - 3] > halfPart) {
        res <- matrix(
                  c(res,
                     length(means), (length(means)) * approxed$k + approxed$b, std.err, approxed$k),
                   ncol = 4, byrow = T)
    } else {
        res <- matrix(res, ncol = 4, byrow = T)
    }
    colnames(res) <- c("x", "y", "std.err", "k")
    res

}

predict <- function(f1, f2,target.probabilty=.95) {
    a1 <- c(f1[1, c("x", "y")])
    a2 <- c(f1[2, c("x", "y")])

    b1 <- c(f2[1, c("x", "y")])
    b2 <- c(f2[2, c("x", "y")])

    if (!are.crossed.right(a1, a2, b1, b2)) {
        res <- c(NA, NA, NA)
    } else {
        a <- sqrt((a1["x"] - b1["x"]) ^ 2 + (a1["y"] - b1["y"]) ^ 2)
        b <- sqrt((a2["x"] - b2["x"]) ^ 2 + (a2["y"] - b2["y"]) ^ 2)
        res.cross.point <- a1 + (a2 - a1) / (1 - (b / a))

        conf.interval.a1 <- a1["std.err"] * 2
        conf.interval.a2 <- a2["std.err"] * 2


        if (a1["y"] > a2["y"]) {
            err.point.a1 <- a1 + c(0, conf.interval.a1)
            err.point.a2 <- a2 - c(0, conf.interval.a2)
        } else {
            err.point.a1 <- a1 - c(0, conf.interval.a1)
            err.point.a2 <- a2 + c(0, conf.interval.a2)
        }

        if (b1["y"] > b2["y"]) {
            err.point.b1 <- b1 + c(0, conf.interval.b1)
            err.point.b2 <- b2 - c(0, conf.interval.b2)
        } else {
            err.point.b1 <- b1 - c(0, conf.interval.a1)
            err.point.b2 <- b2 + c(0, conf.interval.a2)
        }

        err.vect.f1 <- (err.point.a2 - err.point.a1) * (a1["x"] - a2["x"]) / (a1["x"] - res.cross.point["x"])
        err.vect.f2 <- (err.point.b2 - err.point.b1) * (b1["x"] - b2["x"]) / (b1["x"] - res.cross.point["x"])

    }
    res
}


calc.densities <- function(p, d = 60) {
    min.x <- min(p[, "x"])
    max.x <- max(p[, "x"])
    min.y <- min(p[, "y"])
    max.y <- max(p[, "y"])
    step.x <- (max.x - min.x) / d
    step.y <- (max.y - min.y) / d
    cuts.x <- seq(min.x, max.x, by = step.x)
    cuts.y <- seq(min.y, max.y, by = step.y)
    total <- nrow(p)
    a <- cbind(x = p[, "x"],y = p[, "y"], x.c= round((p[, "x"] - min.x) / step.x) * step.x + min.x, y.c = round((p[, "y"] - min.y) / step.y) * step.y + min.y)
    dens <- aggregate(x ~ x.c+ y.c, a, length)      
    colnames(dens) <- c("x", "y", "d")
    dens$d <- dens$d / total * 100
    dens<-dens[dens$d>0.1,]
    m <- matrix(ncol = length(unique(dens$y)), nrow = length(unique(dens$x)))
    rownames(m) <- unique(dens$x)[order(unique(dens$x))]
    colnames(m) <- unique(dens$y)[order(unique(dens$y))]
    apply(dens,MARGIN=1,FUN=function(a) {
        m[as.character(a["x"]),as.character( a["y"])]<<-a["d"]
    })

    list(linear = dens,matr = m)
}

cross.test.check <- function(f1, f2) {

    f1.vect <- f1[, c("x2", "y2")] - f1[, c("x1", "y1")]
    f1.coeff <- cbind(b = f1[, "y1"] - f1.vect[, 2] / f1.vect[, 1] * f1[, "x1"],
                      k = f1.vect[, 2] / f1.vect[, 1])

    f2.vect <- f2[, c("x2", "y2")] - f2[, c("x1", "y1")]
    f2.coeff <- cbind(b = f2[, "y1"] - f2.vect[, 2] / f2.vect[, 1] * f2[, "x1"],
                      k = f2.vect[, 2] / f2.vect[, 1])

    res.cross.point.x <- (f2.coeff[, "b"] - f1.coeff[, "b"]) / (f1.coeff[, "k"] - f2.coeff[, "k"])
    res.cross.point.y <- f1.coeff[, "k"] * res.cross.point.x + f1.coeff[, "b"]
    res.cross.point <- cbind(res.cross.point.x, res.cross.point.y)

    if (nrow(f1) > 1) {
        colnames(res.cross.point) <- c("x", "y")
    } else {
        res.cross.point <- c(x = res.cross.point[1], y = res.cross.point[2])
    }
    res.cross.point
}



make.coeffs <- function(p1, p2) {
    vect <- p2 - p1
    c(b = p1[2] - vect[2] / vect[1] * p1[1], k = vect[2] / vect[1])
}

dens.groups <- function(p) {
    lb <- min(p[!is.na(p[, "d"]), "d"])
    ub <- max(p[!is.na(p[, "d"]), "d"])
    t <- seq(lb, ub, by = (ub - lb) / 6)
    sapply(1:4, FUN =
           function(b) {
               v <- p[p[, "d"] <= t[b] & p[, "d"] > t[b] - (ub - lb) / 4,]

               lines(v[, "x"], v[, "y"],

                    type = "p",
                    col = rgb(1 / b, 1 / b, 1 / b))
           })
}



cross.test <- function(a1.center, a1.sd, a2.center, a2.sd, b1.center, b1.sd, b2.center, b2.sd, test.count = 1000, draw = F, draw.borders = c(0, 10)) {

    f1 <- cbind(
                x1 = seq(a1.center[1], a1.center[1], length = test.count),
                y1 = qnorm(runif(1:test.count, 0, 1), mean = a1.center[2], sd = a1.sd),
                x2 = seq(a2.center[1], a2.center[1], length = test.count),
                y2 = qnorm(runif(1:test.count, 0, 1), mean = a2.center[2], sd = a2.sd))

    f2 <- cbind(
                x1 = seq(b1.center[1], b1.center[1], length = test.count),
                y1 = qnorm(runif(1:test.count, 0, 1), mean = b1.center[2], sd = b1.sd),
                x2 = seq(b2.center[1], b2.center[1], length = test.count),
                y2 = qnorm(runif(1:test.count, 0, 1), mean = b2.center[2], sd = b2.sd))

    t <- cross.test.check(f1, f2)

    if (draw) {
        xs <- c(t[, "x"], a1.center[1], a2.center[1], b1.center[1], b2.center[1])
        xs <- xs[xs > 0 & xs < 10]
        ys <- c(t[, "y"], a1.center[2], a2.center[2], b1.center[2], b2.center[2])
        ys <- ys[ys > 0 & ys < 10]


        f1.vect <- a2.center - a1.center
        f2.vect <- b2.center - b1.center

        plot(0, 0, type = "p",
               xlim = draw.borders,
               ylim = draw.borders, col = "white")
        d <- calc.densities(t)


        lines(t["x"], t["y"], type = "p")
        contour(as.numeric(rownames(d$matr)),as.numeric(colnames(d$matr)), d$matr,add=T,nlevels = 5,col="red")
        abline(coef = make.coeffs(a1.center, a2.center))
        abline(coef = make.coeffs(b1.center, b2.center))

        abline(coef = make.coeffs(a1.center + c(0, a1.sd ), a2.center - c(0, a2.sd)), col = "green", lty = "dotted")
        abline(coef = make.coeffs(a1.center - c(0, a1.sd), a2.center + c(0, a2.sd)), col = "green", lty = "dotted")
        abline(coef = make.coeffs(a1.center + c(0, a1.sd), a2.center + c(0, a2.sd)), col = "red", lty = "dotted")
        abline(coef = make.coeffs(a1.center - c(0, a1.sd), a2.center - c(0, a2.sd)), col = "red", lty = "dotted")

        abline(coef = make.coeffs(b1.center + c(0, b1.sd), b2.center + c(0, b2.sd)), col = "red", lty = "dotted")
        abline(coef = make.coeffs(b1.center - c(0, b1.sd), b2.center + c(0, b2.sd)), col = "green", lty = "dotted")
        abline(coef = make.coeffs(b1.center + c(0, b1.sd), b2.center - c(0, b2.sd )), col = "green", lty = "dotted")
        abline(coef = make.coeffs(b1.center - c(0, b1.sd), b2.center - c(0, b2.sd )), col = "red", lty = "dotted")
    }
    print(density(t))
    print("predicted cross")
    f1 <- matrix(c(a1.center, a2.center, a1.center, a2.center), ncol = 4, byrow = T)
    colnames(f1) <- c("x1", "y1", "x2", "y2")
    f2 <- matrix(c(b1.center, b2.center, b1.center, b2.center), ncol = 4, byrow = T)
    colnames(f2) <- c("x1", "y1", "x2", "y2")
    print(cross.test.check(f1, f2)[1,])

    print("predicted cross2")
    f1 <- matrix(c(a1.center, a2.center, a1.center, a2.center), ncol = 4, byrow = T)

   
    print("crossed.right.percent")
    print(length(t[t[, "x"] > a1.center[1], "x"]) / nrow(t))
    t <- t[t[, "x"] > a1.center[1],]
    print("means")
    print(c(x = mean(t[, "x"]), y = mean(t[, "y"])))
    print("source deviation")
    print(c(a1 = a1.sd, a2 = a2.sd, b1 = b1.sd, b2 = b2.sd))
    print("covariance")
    print(cov(t))
    print("source vectors")
    print(a2.center - a1.center)
    print(b2.center - b1.center)
    print("eigens")
    eig <- eigen(cor(t))
    print(eig)
    cross <- c(mean(t[, "x"]), mean(t[, "y"]))
    axis <- matrix(
        c(cross, cross + eig$vectors[, 1] * eig$values[1] / 2,
        cross, cross - eig$vectors[, 1]  * eig$values[1] / 2,
        cross, cross + eig$vectors[, 2] * eig$values[2] / 2,
        cross, cross - eig$vectors[, 2] * eig$values[2] / 2
        ),
        ncol = 4, byrow = T)
    segments(x0 = axis[,1],y0 = axis[,2],x1=axis[,3],y1 = axis[,4],col="blue")
}



draw.approximation <- function(d, trend = F, err.borders = F) {

    avrg <-
            if (d == 0) {
                cbind(mean = noised, sd = seq(0, 0, length = length(noised)))
            } else {
                ma(noised, d)
            }


    plot(x.values, noised, type = "l", col = "green", ylim = c(min( - avrg[, "sd"] + avrg[, "mean"]), max(avrg[, "sd"] + avrg[, "mean"])))

    if (trend) {
        lines(x = c(trends[, 1], max(x.values)), y = trend.lines[c(trends[, 1], max(x.values))], type = "b", col = "red")
    }

    appr <- swinging.door(avrg[((d - 1) / 2):(nrow(avrg) - (d - 1) / 2),], d)

    lines(x = appr[, "x"] + (d - 1) / 2, y = appr[, "y"], type = "b", col = "blue")
    if (err.borders) {
        lines(x = appr[, "x"] + (d - 1) / 2, y = appr[, "y"] + appr[, "std.err"], type = "l", col = "black", lty = "dashed")
        lines(x = appr[, "x"] + (d - 1) / 2, y = appr[, "y"] - appr[, "std.err"], type = "l", col = "black", lty = "dashed")
    }
}

draw.approximation.std.err <- function(d, trend = F) {

    avrg <-
            if (d == 0) {
                cbind(mean = noised, sd = seq(0, 0, length = length(noised)))
            } else {
                ma(noised, d)
            }
    print(sqrt(mean((avrg[, "mean"] - trend.lines) ^ 2)))
    print(sqrt(mean((avrg[, "sd"]) ^ 2)) / sqrt(d))
    plot(x.values, noised, type = "l", col = "green", ylim = c(min( - avrg[, "sd"] + avrg[, "mean"]), max(avrg[, "sd"] + avrg[, "mean"])))
    if (trend) {
        lines(x = c(trends[, 1], max(x.values)), y = trend.lines[c(trends[, 1], max(x.values))], type = "b", col = "red")
    }

    lines(x.values, avrg[, "mean"], type = "l", col = "black")
    appr <- swinging.door(avrg[((d - 1) / 2):(nrow(avrg) - (d - 1) / 2),], d)
    if (d != 0) {
        lines(x = appr[, "x"] + (d - 1) / 2, 3 * avrg[appr[, "x"], "sd"] / sqrt(d) + appr[, "y"], type = "l", col = "black", lty = "dashed")
        lines(x = appr[, "x"] + (d - 1) / 2, - 3 * avrg[appr[, "x"], "sd"] / sqrt(d) + appr[, "y"], type = "l", col = "black", lty = "dashed")
    }
}



par(mfrow=c(3,1))

#print(sqrt(mean((noised - trend.lines) ^ 2)))
#print(sqrt(mean((ma(noised, 11)[, "sd"]) ^ 2)) / sqrt(11))
#print(sqrt(mean((ma(noised, 101)[, "sd"]) ^ 2)) / sqrt(101))
#print(sqrt(mean((ma(noised, 201)[, "sd"]) ^ 2)) / sqrt(201))
#print(sqrt(mean((ma(noised, 301)[, "sd"]) ^ 2)) / sqrt(301))
#print(sqrt(mean((ma(noised, 501)[, "sd"]) ^ 2)) / sqrt(501))
#print(sqrt(mean((ma(noised, 1001)[, "sd"]) ^ 2)) / sqrt(1001))
#draw.approximation.std.err(5)
#draw.approximation.std.err(51)
draw.approximation.std.err(51)
#draw.approximation(11, err.borders = T, trend = T)
#draw.approximation(51, err.borders = T, trend = T)
#draw.approximation(201, err.borders = T, trend = T)


#par(mfrow = c(1, 1))
#cross.result<-cross.test(c(0,0),1,c(15,15),.2,c(0,15),1,c(20,0),1,test.count = 10000000,draw = T,draw.borders = c(5,10))


