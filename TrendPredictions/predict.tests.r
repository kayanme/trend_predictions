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
    a <- cbind(x = p[, "x"], y = p[, "y"], x.c = round((p[, "x"] - min.x) / step.x) * step.x + min.x, y.c = round((p[, "y"] - min.y) / step.y) * step.y + min.y)
    dens <- aggregate(x ~ x.c + y.c, a, length)
    colnames(dens) <- c("x", "y", "d")
    dens$d <- dens$d / total * 100
    dens <- dens[dens$d > 0.1,]
    m <-
    outer(unique(dens$x), unique(dens$y), FUN = function(xl, yl) {
        unlist(apply(cbind(xl, yl), MARGIN = 1, FUN = function(xy) {
            a <- dens[dens$x == xy[1] & dens$y == xy[2],]$d
            if (length(a) > 0) {
                c(a)
            } else {
                0
            }
        }))

    })

    rownames(m) <- unique(dens$x)[order(unique(dens$x))]
    colnames(m) <- unique(dens$y)[order(unique(dens$y))]


    list(linear = dens, matr = m)
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

prepare.checks <- function(a1.center, a1.sd, a2.center, a2.sd, b1.center, b1.sd, b2.center, b2.sd, test.count = 1000) {
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

    cross.test.check(f1, f2)
}

cross.test.data <- function(test.data) {
    print("means")
    print(c(x = mean(test.data[, "x"]), y = mean(test.data[, "y"])))
    print("covariance")
    print(cov(test.data))
    print("eigens")
    eig <- eigen(cov(test.data))
    print(eig)
}

cross.test.real <- function(a1.center, a1.sd, a2.center, a2.sd, b1.center, b1.sd, b2.center, b2.sd, test.data) {
    to.draw <- test.data
    f1.vect <- a2.center - a1.center
    f2.vect <- b2.center - b1.center

  
    lines(to.draw[, "x"], to.draw[, "y"], type = "p", col = "black")   
    abline(coef = make.coeffs(a1.center, a2.center))
    abline(coef = make.coeffs(b1.center, b2.center))

    abline(coef = make.coeffs(a1.center + c(0, a1.sd), a2.center - c(0, a2.sd)), col = "green", lty = "dotted")
    abline(coef = make.coeffs(a1.center - c(0, a1.sd), a2.center + c(0, a2.sd)), col = "green", lty = "dotted")
    abline(coef = make.coeffs(a1.center + c(0, a1.sd), a2.center + c(0, a2.sd)), col = "red", lty = "dotted")
    abline(coef = make.coeffs(a1.center - c(0, a1.sd), a2.center - c(0, a2.sd)), col = "red", lty = "dotted")

    abline(coef = make.coeffs(b1.center + c(0, b1.sd), b2.center + c(0, b2.sd)), col = "red", lty = "dotted")
    abline(coef = make.coeffs(b1.center - c(0, b1.sd), b2.center + c(0, b2.sd)), col = "green", lty = "dotted")
    abline(coef = make.coeffs(b1.center + c(0, b1.sd), b2.center - c(0, b2.sd)), col = "green", lty = "dotted")
    abline(coef = make.coeffs(b1.center - c(0, b1.sd), b2.center - c(0, b2.sd)), col = "red", lty = "dotted")
}

cross.test.real.density <- function(test.data) {
    d <- calc.densities(to.draw, 100)
    contour(as.numeric(rownames(d$matr)), as.numeric(colnames(d$matr)), d$matr, add = T, nlevels = 5, col = "blue")
}

cross.test.real.axis <- function(test.data) {
    cross <- c(mean(test.data[, "x"]), mean(test.data[, "y"]))
    axis <- matrix(
    c(cross, cross + eig$vectors[, 1] * eig$values[1] * 10,
    cross, cross - eig$vectors[, 1] * eig$values[1] * 10,
    cross, cross + eig$vectors[, 2] * eig$values[2] * 10,
    cross, cross - eig$vectors[, 2] * eig$values[2] * 10),
    ncol = 4, byrow = T)
    segments(x0 = axis[, 1], y0 = axis[, 2], x1 = axis[, 3], y1 = axis[, 4], col = "blue")
}

shift <- function(x1, y1, x2, y2, x3, y3, x4, y4) {
    J <- cbind(
      c(y1 - y2, x3 - x4),
      c(x2 - x1, x4 - x3))
    J <- det(J)
    det1 <- y1 * x2 - y2 * x1
    det2 <- y3 * x4 - y4 * x3

    function(x = NULL, y = NULL, point = NULL) {
        if (is.null(point)) {
            shifted.x <- (x3 - x4) * x - (x1 - x2) * y + (x3 - x4) *  det1 - (x1 - x2) *  det2
            shifted.x <- shifted.x / J

            shifted.y <- (y3 - y4) * x - (y1 - y2) * y + (y3 - y4) * det1 - (y1 - y2) * det2
            shifted.y <- shifted.y / J

            cbind(x = shifted.x, y = shifted.y)
        } else {
            c((x3 - x4) * (point[1] + det1)
            - (x1 - x2) * (point[2] + det2),
                  (y3 - y4) * (point[1] + det1)
                  - (y1 - y2) * (point[2] + det2)) / J
        }
    }
}

cross.test.data.shifted <- function(test.data) {
    print("means")
    print(c(x = mean(test.data[, "x"]), y = mean(test.data[, "y"])))
    cfs <- make.coeffs(c(test.data[1,]), c(test.data[100,]))
    print(cfs)
}

cross.test.shifted <- function(a1.center, a1.sd, a2.center, a2.sd, b1.center, b1.sd, b2.center, b2.sd, test.data) {

    to.draw <- test.data

    shift.fun <- shift(a1.center[1], a1.center[2], a2.center[1], a2.center[2],
                        b1.center[1], b1.center[2], b2.center[1], b2.center[2])

    f1.vect <- a2.center - a1.center
    f2.vect <- b2.center - b1.center

    shifted <- shift.fun(to.draw[, "x"], to.draw[, "y"])

    lines(shifted[, 1], shifted[, 2], type = "p", col = "black")


    lines(mean(shifted[, 1]), mean(shifted[, 2]), type = "p", col = "red")
    draw.shifted <- function(p1, p2, ...) {
        s <- shift.fun(point = p1)
        s2 <- shift.fun(point = p2)
        if (s[1] == s2[1]) {
            abline(v = s[1], ...)
        } else {
            abline(coef = make.coeffs(s, s2), ...)
        }


    }
    cross.test.data.shifted(shifted)

    draw.shifted(a1.center, a2.center)
    draw.shifted(b1.center, b2.center)


    draw.shifted(a1.center + c(0, a1.sd), a2.center - c(0, a2.sd), col = "green", lty = "dotted")
    draw.shifted(a1.center - c(0, a1.sd), a2.center + c(0, a2.sd), col = "green", lty = "dotted")
    draw.shifted(a1.center + c(0, a1.sd), a2.center + c(0, a2.sd), col = "red", lty = "dotted")
    draw.shifted(a1.center - c(0, a1.sd), a2.center - c(0, a2.sd), col = "red", lty = "dotted")

    draw.shifted(b1.center + c(0, b1.sd), b2.center + c(0, b2.sd), col = "red", lty = "dotted")
    draw.shifted(b1.center - c(0, b1.sd), b2.center + c(0, b2.sd), col = "green", lty = "dotted")
    draw.shifted(b1.center + c(0, b1.sd), b2.center - c(0, b2.sd), col = "green", lty = "dotted")
    draw.shifted(b1.center - c(0, b1.sd), b2.center - c(0, b2.sd), col = "red", lty = "dotted")



}