#    Copyright (C) 2016 University of Southern California and
#             Chao Deng and Andrew D. Smith and Timothy Daley
#
#    Authors: Chao Deng
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

### initial settings of two parameters size and mu in a negative binomial
### distribution for a numeric optimal searching function optim in R
SIZE.INIT <- 1
MU.INIT <- 0.5

### termination conditions for EM algorithm
TOLERANCE <- 1e-10
ITER.TOLERANCE <- 1e5


### density function of a truncated zero negative binomial distribution
### size and mu are two parameters for the negative binomial
zerotruncated.dnbinom <- function(x, size, mu, log = FALSE)
{
    ## the density of x in negative binomial
    p <- dnbinom(x, size = size, mu = mu, log = log)
    
    ## set zeros in x with zero probability
    if (log == FALSE) {
        p[ which(x == 0) ] <- 0
    } else {
        p[ which(x == 0) ] <- -Inf
    }
    
    ## the density of non-zero in negative binomial
    q <- 1 - dnbinom(0, size = size, mu = mu)
    
    ## normalize all non-zero values in negrative binomial to generate ZTNB
    if (log == FALSE) {
        return( p/q )
    } else {
        return( p - log(q) )
    }
}


### zerotruncated negative loglikelihood
zerotruncated.minus.log.likelihood <- function(n, size, mu)
{
    prob <- zerotruncated.dnbinom(n[, 1], size, mu, log = TRUE)
    
    ## negative loglikelihood
    prob <- -prob
    return( prob %*% n[, 2] )
}


### calculate the negative binomial loglikelihood
### zero.items is number of items unobserved
### size and mu are parameters in a negative binomial distribution
nb.loglikelihood <- function(n, zero.items, size, mu)
{
    ## likelihood of nonzero terms
    log.prob <- dnbinom(n[, 1], size = size, mu = mu, log = TRUE)
    loglikelihood <- log.prob %*% n[, 2]
    
    ## add items with zero count
    log.zero.prob <- dnbinom(0, size = size, mu = mu, log = TRUE)
    loglikelihood <- loglikelihood + zero.items * log.zero.prob
    
    return(loglikelihood)
}


### EM algorithm to fit the histogram with a negative binomial distribution
### hist only includes information for observation
### the number of unobserved items is missing data
preseqR.ztnb.em <- function(n, size=SIZE.INIT, mu=MU.INIT)
{
    checking.hist(n)
    
    n[, 2] <- as.numeric(n[, 2])
    ## setting the number of unobserved items as 0
    zero.prob <- exp(dnbinom(0, size = size, mu = mu, log = TRUE))
    
    ## estimate the total number of distinct items
    observed.items <- sum(n[, 2])
    L <- observed.items/( 1 - zero.prob )
    
    ## expected the number of unobservations
    zero.items <- L*zero.prob
    
    ## estimated mean and variance
    m <- (n[, 1] %*% n[, 2]) / L
    v <- ( (n[, 1] - m)^2 %*% n[, 2] + m^2 * zero.items )/(L - 1)
    
    ## target function f
    f <- function(x) {
        return( -nb.loglikelihood(n, zero.items, size = x, mu = m)/L )
    }
    
    ## derivative of f
    gr <- function(x)
    {
        first.term <- ( digamma(x) * zero.items +
                            digamma(n[, 1] + x) %*% n[, 2] )/L
        second.term <- digamma(x)
        third.term <- log(x) - log(x + m)
        result <- first.term - second.term + third.term
        # f is negative loglikelihood
        return(-result)
    }
    
    ## estimate size and mu based on first and second moments
    if (v > m) {
        res <- optim(m^2 / (v - m), f, gr, method = "L-BFGS-B",
                     lower = 0.0001, upper = 10000)
    } else {
        res <- optim(size, f, gr, method = "L-BFGS-B",
                     lower = 0.0001, upper = 10000)
    }
    
    ## count the times of iteration
    iter <- as.double(1)
    
    ## initialize the negative loglikelihood
    loglikelihood.pre <- Inf
    
    ## zerotruncated loglikelihood
    loglikelihood <- zerotruncated.minus.log.likelihood(n, res$par, m)
    
    ## EM algorithm
    while (( loglikelihood.pre - loglikelihood )/observed.items > TOLERANCE &&
           iter < ITER.TOLERANCE)
    {
        ## update negative loglikelihood
        loglikelihood.pre <- loglikelihood
        
        ## update parameters
        size <- res$par
        mu <- m
        
        ### E-step: estimate the number of unobserved items
        
        ## update the probility an item unobserved
        zero.prob <- exp(dnbinom(0, size = size, mu = mu, log = TRUE))
        
        ## estimate the total number of distinct items
        L <- observed.items/( 1 - zero.prob )
        
        ## update expected number of unobserved items
        zero.items <- L*zero.prob
        
        ## estimated mean and variance
        m <- (n[, 1] %*% n[, 2])/L
        v <- ( (n[, 1] - m)^2 %*% n[, 2] + m^2 * zero.items )/(L - 1)
        
        ### M step: estimate the parameters size and mu
        if (v > m) {
            res <- optim(m^2 / (v - m), f, gr, method = "L-BFGS-B",
                         lower = 0.0001, upper = 10000)
        } else {
            res <- optim(size, f, gr, method = "L-BFGS-B",
                         lower = 0.0001, upper = 10000)
        }
        iter <- iter + 1
        ## zerotruncated loglikelihood
        loglikelihood <- zerotruncated.minus.log.likelihood(n, res$par, m)
    }
    return(list(size = size, mu = mu, loglik = -loglikelihood.pre))
}


## fitting the negative binoimal distribution to the data by EM algorithm
## r is a vector of frequencies
## return an estimator by ZTNB
ztnb.mincount <- function(n, r=1, size=SIZE.INIT, mu=MU.INIT)
{
    checking.hist(n)
    
    n[, 2] <- as.numeric(n[, 2])
    total.sample <- n[, 1] %*% n[, 2]
    distinct <- sum(n[, 2])
    
    ## estimate parameters
    opt <- preseqR.ztnb.em(n, size, mu)
    size <- opt$size
    mu <- opt$mu
    
    ## the probability of being sampled in the initial experiment
    p <- 1 - dnbinom(0, size = size, mu = mu)
    
    ## L is the estimated total number of distinct items
    L <- distinct/p
    
    f.mincount <- function(t) {
        L * pnbinom(r - 1, size=size, mu=mu*t, lower.tail=FALSE)
    }
    f.mincount(1); f.mincount
}





#    Copyright (C) 2016 University of Southern California and
#             Chao Deng and Andrew D. Smith and Timothy Daley
#
#    Authors: Chao Deng
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

## Two roots are the same if the difference is less than the PRECISION
PRECISION <- 1e-3

### interpolating for species accumulation curve with minimum count
### ss step size
### n two-column histogram
### r minimum count
preseqR.interpolate.mincount <- function(ss, n, r=1)
{
    checking.hist(n)
    
    n[, 2] <- as.numeric(n[, 2])
    ## total individuals captured
    total.sample <- n[, 1] %*% n[, 2]
    N <- total.sample
    
    ## total species
    initial.distinct <- sum(n[, 2])
    step.size <- as.double(ss)
    
    ## l is the number of sampled points for interpolation
    l <- N / step.size
    
    ## if the sample size is larger than the size of experiment or 
    ## the step size is too small, return NULL
    if (l < 1 || ss < 1 || r < 1)
        return()
    ## if the sample size is the size of the experiment
    ## count the number of species observed r or more times
    else if (l == 1) {
        index <- which(n[, 1] >= r)
        result <- matrix(c(step.size, sum(n[index, 2])), ncol = 2, byrow = FALSE)
        colnames(result) <- c('sample.size', 'interpolation')
        return(result)
    }
    
    ## explicitly calculating the expected species observed at least r times
    ## based on sampling without replacement
    ## see K.L Heck 1975
    ## N total individuals
    ## S the number of species
    ## size the size of the subsample
    expect.distinct <- function(n, N, size, S, r) {
        denom <- lchoose(N, size)
        p <- sapply(n[, 1], function(x) {
            sum(exp(lchoose(N - x, size - 0:(r-1)) + lchoose(x, 0:(r-1)) - denom))})
        return(S - p %*% n[, 2])
    }
    
    ## sample sizes
    x <- step.size * ( 1:l )
    
    ## calculate the number of distinct reads based on each sample size
    yield.estimates <- sapply(x, function(x) {
        expect.distinct(n, N, x, initial.distinct, r)})
    
    ## put size and yield together into a matrix
    result <- matrix(c(x, yield.estimates), ncol = 2, byrow = FALSE)
    colnames(result) <- c('sample.size', 'interpolation')
    
    return(result)
}


### power series based on count frequencies starting from frequency j
### when j = 1, it is the power series expansion of E(S_1(t)) / t at t = 1
### the maximum number of terms
generating.ps <- function(n, mt, j=1) {
    if (j >= max(n[, 1])) return(NULL)
    ## transform a histogram into a vector of frequencies
    hist.count <- vector(length=max(n[, 1]), mode="numeric")
    hist.count[n[, 1]] <- n[, 2]
    
    ## shift to required count frequencies
    hist.count <- hist.count[j: length(hist.count)]
    
    PS.coeffs <- sum(hist.count)
    change.sign <- 0
    
    ## preserve extra precision mt+1
    for (i in 1:(min(mt+1, length(hist.count)))) {
        PS.coeffs <- c(PS.coeffs, 
                       (-1)^change.sign * hist.count[i] - PS.coeffs[length(PS.coeffs)])
        change.sign <- change.sign + 1
    }
    
    ## truncate at coefficients where it is zero
    zero.index <- which(PS.coeffs == 0)
    if (length(zero.index) > 0) {
        PS.coeffs[1:(min(zero.index) - 1)]
    } else {
        PS.coeffs
    }
}


## species accum curves based on parital fraction expansion
## the function is used whenever count frequency 1 is unavaible or the sample
## size is saturated
## using count frequencies starting from a given count frequency instead of 1
## when start.freq = 1, it is identical to the function preseqR.pf.mincount
## CHAO: save for a rainy day
general.ds.mincount <- function(n, r=1, mt=20, start.freq=1)
{
    # check the input format of the histogram
    checking.hist(n)
    
    n[, 2] <- as.numeric(n[, 2])
    ## constructing the power series
    PS.coeffs <- generating.ps(n, j=start.freq, mt=mt)
    
    if (is.null(PS.coeffs)) {
        write("the size of the initial experiment is insufficient", stderr())
        return(NULL)
    }
    
    ## constrain the continued fraction approximation with even degree 
    ## asymptotically ~ C / t
    mt <- min(mt, length(PS.coeffs))
    PS.coeffs <- PS.coeffs[ 1:mt ]
    
    ## check whether sample size is sufficient
    if (mt < 2)
    {
        m <- paste("max count before zero is les than min required count (2)",
                   " sample not sufficiently deep or duplicates removed", sep = ',')
        write(m, stderr())
        return(NULL)
    }
    
    ## construct the continued fraction approximation to the power seies
    cf <- ps2cfa(coef=PS.coeffs, mt=mt)
    rf <- cfa2rf(CF=cf)
    ## the length of cf could be less than mt
    ## even if ps do not have zero terms, coefficients of cf may have
    mt <- length(cf)
    ## select rational function approximants [M-1/M] m=2M
    ## asymptotically ~ C / t
    mt <- mt - (mt %% 2)
    valid.estimator <- FALSE
    m <- mt
    while (valid.estimator == FALSE) {
        
        rfa <- rf2rfa(RF=rf, m=m)
        ## solving roots
        numer.roots <- solve(rfa[[1]])
        denom.roots <- solve(rfa[[2]])
        ## seperating roots by their real parts
        numer.roots.neg <- numer.roots[which(Re(numer.roots) < 0)]
        numer.roots.pos <- numer.roots[which(Re(numer.roots) >= 0)]
        denom.roots.neg <- denom.roots[which(Re(denom.roots) < 0)]
        denom.roots.pos <- denom.roots[which(Re(denom.roots) >= 0)]
        
        ## record roots in the numerator that are significantly similar to
        ## roots in the denominator
        tmp.roots <- c()
        
        ## simplify the rational function approximation
        ## two roots are same if the difference is less than the 
        ## predefined PRECISION
        if (length(numer.roots.pos) > 0) {
            for (i in 1:length(numer.roots.pos)) {
                if (length(denom.roots.pos) > 0) {
                    d <- Mod(denom.roots.pos - numer.roots.pos[i])
                    for (j in 1:length(d)) {
                        if (d[j] < PRECISION) {
                            denom.roots.pos <- denom.roots.pos[-j]
                            tmp.roots <- c(tmp.roots, numer.roots.pos[i])
                            break
                        }
                    }
                }
            }
        }
        
        ## roots in simplified RFA
        numer.roots <- numer.roots[!numer.roots %in% tmp.roots]
        denom.roots <- c(denom.roots.neg, denom.roots.pos)
        
        ## convert roots from t - 1 to t
        roots <- denom.roots + 1
        ## pacman rule checking
        if (length(which(roots == 0)) || length(which(Re(roots) > 0))) {
            m <- m - 2
            next
        } else {
            poly.numer <- as.function(poly.from.roots(numer.roots))
            l <- length(denom.roots)
            ## treat polynomials in the rational function to be monic
            ## the difference to the original RFA is a multiplier C
            
            ## c_i in the estimator
            coef <- sapply(1:l, function(x) {
                poly.numer(denom.roots[x]) / prod(denom.roots[x] - denom.roots[-x])})
            ## calculate the constant C
            C <- coef(rfa[[1]])[length(coef(rfa[[1]]))] / 
                coef(rfa[[2]])[length(coef(rfa[[2]]))]
            ## species accum curves with minimum count r
            ## using parital fraction expansion
            denom.roots <- denom.roots + 1
            coef <- coef * C
            ## modify the coefficients
            coef <- coef * (1 - denom.roots)^(start.freq - 1)
            ## check whether the estimator is non-decreased                             
            deriv.f <- function(t) {
                Re(sapply(t, function(x) {-(coef*denom.roots) %*% ( 1 / ((x-denom.roots)^2))}))} 
            if (length(which( deriv.f(seq(0.05, 100, by=0.05)) < 0 ) != 0)) {
                m <- m - 2
                next
            }
            f.mincount <- function(t) {
                sapply(r, function(x) {
                    Re(coef %*% (t / (t - denom.roots))^x)})}
            f.mincount(1)
            valid.estimator <- TRUE
        }
    }
    ## remove M, M.adjust in the future
    list(FUN=f.mincount, M=m / 2, M.adjust=length(denom.roots), FUN.elements=list(coef=coef, roots=denom.roots))
}

### power series based on count frequencies starting from frequency j
### when j = 1, it is the power series expansion of E(S_1(t)) / t at t = 1
### the maximum number of terms
generating.ps <- function(n, mt, j=1) {
    if (j >= max(n[, 1])) return(NULL)
    ## transform a histogram into a vector of frequencies
    hist.count <- vector(length=max(n[, 1]), mode="numeric")
    hist.count[n[, 1]] <- n[, 2]
    
    ## shift to required count frequencies
    hist.count <- hist.count[j: length(hist.count)]
    
    PS.coeffs <- sum(hist.count)
    change.sign <- 0
    
    ## preserve extra precision mt+1
    for (i in 1:(min(mt+1, length(hist.count)))) {
        PS.coeffs <- c(PS.coeffs, 
                       (-1)^change.sign * hist.count[i] - PS.coeffs[length(PS.coeffs)])
        change.sign <- change.sign + 1
    }
    
    ## truncate at coefficients where it is zero
    zero.index <- which(PS.coeffs == 0)
    if (length(zero.index) > 0) {
        PS.coeffs[1:(min(zero.index) - 1)]
    } else {
        PS.coeffs
    }
}


## species accum curves based on parital fraction expansion
## the function is used whenever count frequency 1 is unavaible or the sample
## size is saturated
## using count frequencies starting from a given count frequency instead of 1
## when start.freq = 1, it is identical to the function preseqR.pf.mincount
## CHAO: save for a rainy day
general.ds.mincount <- function(n, r=1, mt=20, start.freq=1)
{
    # check the input format of the histogram
    checking.hist(n)
    
    n[, 2] <- as.numeric(n[, 2])
    ## constructing the power series
    PS.coeffs <- generating.ps(n, j=start.freq, mt=mt)
    
    if (is.null(PS.coeffs)) {
        write("the size of the initial experiment is insufficient", stderr())
        return(NULL)
    }
    
    ## constrain the continued fraction approximation with even degree 
    ## asymptotically ~ C / t
    mt <- min(mt, length(PS.coeffs))
    PS.coeffs <- PS.coeffs[ 1:mt ]
    
    ## check whether sample size is sufficient
    if (mt < 2)
    {
        m <- paste("max count before zero is les than min required count (2)",
                   " sample not sufficiently deep or duplicates removed", sep = ',')
        write(m, stderr())
        return(NULL)
    }
    
    ## construct the continued fraction approximation to the power seies
    cf <- ps2cfa(coef=PS.coeffs, mt=mt)
    rf <- cfa2rf(CF=cf)
    ## the length of cf could be less than mt
    ## even if ps do not have zero terms, coefficients of cf may have
    mt <- length(cf)
    ## select rational function approximants [M-1/M] m=2M
    ## asymptotically ~ C / t
    mt <- mt - (mt %% 2)
    valid.estimator <- FALSE
    m <- mt
    while (valid.estimator == FALSE) {
        
        rfa <- rf2rfa(RF=rf, m=m)
        ## solving roots
        numer.roots <- solve(rfa[[1]])
        denom.roots <- solve(rfa[[2]])
        ## seperating roots by their real parts
        numer.roots.neg <- numer.roots[which(Re(numer.roots) < 0)]
        numer.roots.pos <- numer.roots[which(Re(numer.roots) >= 0)]
        denom.roots.neg <- denom.roots[which(Re(denom.roots) < 0)]
        denom.roots.pos <- denom.roots[which(Re(denom.roots) >= 0)]
        
        ## record roots in the numerator that are significantly similar to
        ## roots in the denominator
        tmp.roots <- c()
        
        ## simplify the rational function approximation
        ## two roots are same if the difference is less than the 
        ## predefined PRECISION
        if (length(numer.roots.pos) > 0) {
            for (i in 1:length(numer.roots.pos)) {
                if (length(denom.roots.pos) > 0) {
                    d <- Mod(denom.roots.pos - numer.roots.pos[i])
                    for (j in 1:length(d)) {
                        if (d[j] < PRECISION) {
                            denom.roots.pos <- denom.roots.pos[-j]
                            tmp.roots <- c(tmp.roots, numer.roots.pos[i])
                            break
                        }
                    }
                }
            }
        }
        
        ## roots in simplified RFA
        numer.roots <- numer.roots[!numer.roots %in% tmp.roots]
        denom.roots <- c(denom.roots.neg, denom.roots.pos)
        
        ## convert roots from t - 1 to t
        roots <- denom.roots + 1
        ## pacman rule checking
        if (length(which(roots == 0)) || length(which(Re(roots) > 0))) {
            m <- m - 2
            next
        } else {
            poly.numer <- as.function(poly.from.roots(numer.roots))
            l <- length(denom.roots)
            ## treat polynomials in the rational function to be monic
            ## the difference to the original RFA is a multiplier C
            
            ## c_i in the estimator
            coef <- sapply(1:l, function(x) {
                poly.numer(denom.roots[x]) / prod(denom.roots[x] - denom.roots[-x])})
            ## calculate the constant C
            C <- coef(rfa[[1]])[length(coef(rfa[[1]]))] / 
                coef(rfa[[2]])[length(coef(rfa[[2]]))]
            ## species accum curves with minimum count r
            ## using parital fraction expansion
            denom.roots <- denom.roots + 1
            coef <- coef * C
            ## modify the coefficients
            coef <- coef * (1 - denom.roots)^(start.freq - 1)
            ## check whether the estimator is non-decreased                             
            deriv.f <- function(t) {
                Re(sapply(t, function(x) {-(coef*denom.roots) %*% ( 1 / ((x-denom.roots)^2))}))} 
            if (length(which( deriv.f(seq(0.05, 100, by=0.05)) < 0 ) != 0)) {
                m <- m - 2
                next
            }
            f.mincount <- function(t) {
                sapply(r, function(x) {
                    Re(coef %*% (t / (t - denom.roots))^x)})}
            f.mincount(1)
            valid.estimator <- TRUE
        }
    }
    ## remove M, M.adjust in the future
    list(FUN=f.mincount, M=m / 2, M.adjust=length(denom.roots), FUN.elements=list(coef=coef, roots=denom.roots))
}


## nonparametric approach Deng & Smith 2016
ds.mincount.bootstrap <- function(n, r=1, mt=20, times=100)
{
    n[, 2] <- as.numeric(n[, 2])
    ## total individuals
    total <- n[, 1] %*% n[, 2]
    
    ## returned function
    f.mincount <- vector(length=times, mode="list")
    
    ds.estimator <- function(n, r, mt, t.scale) {
        f <- ds.mincount(n, r=r, mt=mt)
        if (f$M == 1) {
            f <- ztnb.mincount(n, r=r)
            function(t) {f(t * t.scale)}
        } else {
            function(t) {f$FUN(t * t.scale)}
        }
    }
    
    while (times > 0) {
        n.bootstrap <- matrix(c(n[, 1], rmultinom(1, sum(n[, 2]), n[, 2])), ncol=2)
        total.bootstrap <- n.bootstrap[, 1] %*% n.bootstrap[, 2]
        t.scale <- total / total.bootstrap
        f <-  ds.estimator(n.bootstrap, r=r, mt=mt, t.scale=t.scale) 
        
        f.mincount[[times]] <- f
        ## prevent later binding!!!
        f.mincount[[times]](1)
        times <- times - 1
    }
    f.estimator <- ds.mincount(n=n, r=r, mt=mt)
    if (length(r) == 1) {
        median.estimators <- function(t) {median( sapply(f.mincount, function(x) x(t)) )}
        var.estimator <- function(t) {var( sapply(f.mincount, function(x) x(t)) )}
    } else {
        median.estimators <- function(t) {apply(sapply(f.mincount, function(x) x(t)), FUN=median, MARGIN=1)}
        var.estimator <- function(t) {apply(sapply(f.mincount, function(x) x(t)), FUN=var, MARGIN=1)}
    }
    ## prevent later binding!!!
    f.estimator$FUN(1); median.estimators(1); var.estimator(1)
    return(list(FUN.nobootstrap=f.estimator, FUN.bootstrap=median.estimators, var=var.estimator))
}

ds.mincount <- function(n, r=1, mt=20)
{
    checking.hist(n)
    
    n[, 2] <- as.numeric(n[, 2])
    
    ## constructing the power series
    PS.coeffs <- generating.ps(n, mt=mt, j=1)
    
    if (is.null(PS.coeffs)) {
        write("the size of the initial experiment is insufficient", stderr())
        return(NULL)
    }
    
    ## use only the first mt terms in the power series
    mt <- min(mt, length(PS.coeffs))
    PS.coeffs <- PS.coeffs[ 1:mt ]
    
    ## check whether sample size is sufficient
    if (mt < 2)
    {
        m <- paste("max count before zero is less than min required count (2)",
                   " sample not sufficiently deep or duplicates removed", sep = ',')
        write(m, stderr())
        return(NULL)
    }
    
    ## construct the continued fraction approximation to the power seies
    cf <- ps2cfa(coef=PS.coeffs, mt=mt)
    rf <- cfa2rf(CF=cf)
    ## the length of cf could be less than mt
    ## even if ps do not have zero terms, coefficients of cf may have
    mt <- length(cf)
    ## select rational function approximants [M-1/M] m=2M
    ## asymptotically ~ C / t
    mt <- mt - (mt %% 2)
    valid.estimator <- FALSE
    m <- mt
    while (valid.estimator == FALSE && m >= 2) {
        
        rfa <- rf2rfa(RF=rf, m=m)
        ## solving roots
        numer.roots <- solve(rfa[[1]])
        denom.roots <- solve(rfa[[2]])
        
        ## finite
        if (any(!is.finite(c(numer.roots, denom.roots)))) {
            m = m - 2
            next;
        }
        
        ## record roots in the numerator that are significantly similar to
        ## roots in the denominator
        tmp.roots <- c()
        
        ## simplify the rational function approximation
        ## two roots are same if the difference is less than the 
        ## predefined PRECISION
        if (length(denom.roots) > 0) {
            for (i in 1:length(denom.roots)) {
                if (length(numer.roots) > 0) {
                    d <- Mod(denom.roots[i] - numer.roots)
                    ind <- which.min(d)
                    if (d[ind] < PRECISION) {
                        numer.roots <- numer.roots[-ind]
                        tmp.roots <- c(tmp.roots, denom.roots[i])
                    }
                }
            }
        }
        
        ## roots in simplified RFA
        denom.roots <- denom.roots[!denom.roots %in% tmp.roots]
        
        ## convert roots from t - 1 to t
        roots <- denom.roots + 1
        
        ## pacman rule checking
        if (any(Re(roots) >= 0)) {
            m <- m - 2
            next
        } else {
            if (length(numer.roots) == 0) {
                poly.numer <- as.function(polynomial(1))
            } else {
                poly.numer <- as.function(poly.from.roots(numer.roots))
            }
            l <- length(denom.roots)
            ## treat polynomials in the rational function to be monic
            ## the difference to the original RFA is a multiplier C
            
            ## c_i in the estimator
            coefs <- sapply(1:l, function(x) {
                poly.numer(denom.roots[x]) / prod(denom.roots[x] - denom.roots[-x])})
            ## calculate the constant C
            C <- coef(rfa[[1]])[length(coef(rfa[[1]]))] / 
                coef(rfa[[2]])[length(coef(rfa[[2]]))]
            coefs <- coefs * C
            
            ## check whether the estimator is non-decreased
            ## NOTE: it only checks for t >= 1 !!!
            
            deriv.f <- function(t) {
                Re(sapply(t, function(x) {-(coefs*roots) %*% ( 1 / ((x-roots)^2))}))}
            if (any( deriv.f(seq(1, 100, by=0.05)) < 0 )) {
                m <- m - 2
                next
            } else {
                f.mincount <- function(t) {
                    sapply(r, function(x) {
                        Re(coefs %*% (t / (t - roots))^x)})}
                f.mincount(1)
                valid.estimator <- TRUE
            }
        }
    }
    ## remove M, M.adjust in the future)
    if (valid.estimator == TRUE) {
        return(list(FUN=f.mincount, M=m / 2, M.adjust=length(roots), FUN.elements=list(coefs=coefs, roots=roots)))
    } else {
        ## the case m = 0
        f.mincount <- function(t) {
            sapply(r, function(x) sum(n[, 2]))}
        return(list(FUN=f.mincount, M=1, M.adjust=1, FUN.elements=list(coefs=sum(n[, 2]), roots=0)))
    }
}


#    Copyright (C) 2016 University of Southern California and
#             Chao Deng and Andrew D. Smith and Timothy Daley
#
#    Authors: Chao Deng
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


### checking the input histogram in an appropariat format
checking.hist <- function(n)
{
    if (ncol(n)!=2 || is.numeric(n[,1])==FALSE || is.numeric(n[,2])==FALSE) {
        stop("Input must be a two-column matrix")
    }
    ## the first column is the frequencies of observed items
    freq <- n[, 1]
    
    ## the second column is the number of observed distinct items for each
    ## frequency
    number.items <- n[, 2]
    
    ## check whether frequencies are at least one and the histogram is sorted
    for (i in 1:length(freq))
        if (freq[i] <= 0 || freq[i] != floor(freq[i])) {
            stop("The first column must be positive integers!")
        } else if (number.items[i] < 0) {
            stop("The second column must be non negative")
        }
    else {
        if (i > 1 && freq[i - 1] >= freq[i])
            stop("The first column is not sorted in the ascending order")
    }
    
    return(n)
}

### check determinants of matrix M_{m-1,m-1},M_{m-1,m},M_{m,m-1},M_{m,m}
checking.matrix.det <- function(n, m) {
    ps <- generating.ps(n, j=1, mt=2*m + 1)
    ps <- c(0, ps)
    matrix.dets <- vector(length=4, mode="numeric")
    count <- 1
    for (i in (m-1):m)
        for (j in (m-1):m) {
            pade.matrix <- sapply(1:(i+1), function(x) {
                start <- j + 1 + x
                end <- j - i + 1 + x
                indexes <- seq(start, end, -1)
                ps[indexes]})
            
            matrix.dets[count] <- det(pade.matrix / max(abs(pade.matrix)))
            count <- count + 1
        }
    matrix.dets
}



## sampling without replacement
## n frequencies counts
nonreplace.sampling <- function(size, n)
{
    ## make sure frequencies are integers
    n[, 2] <- floor(n[, 2])
    ## the number of distinct items
    distinct <- sum(n[, 2])
    
    ## identifier for each distinct item
    ind <- 1:distinct
    
    ## the size of each read in the library
    N <- rep(n[, 1], n[, 2])
    
    ## construct a sample space X 
    ## the whole library represents by its indexes. If a read presents t
    ## times in the library, its indexes presents t times in X
    X <- rep(ind, N)
    
    return(sample(X, size, replace = FALSE))
}


## sampling without replacement
## input frequencies counts; output subsample as a frequencies counts
preseqR.nonreplace.sampling <- function(size, n)
{
    ## check the input histogram file
    checking.hist(n)
    ## sub sampling
    X <- nonreplace.sampling(size, n)
    ## record the freq of each sampled species
    freq.counts <- hist(X, breaks=0:max(X), plot=FALSE)$count
    ## frequencies counts; frequency 0 excluded
    T <- hist(freq.counts, breaks=-1:max(freq.counts), plot=FALSE)$counts[-1]
    matrix(c(which(T != 0), T[which(T != 0)]), byrow = FALSE, ncol=2)
}


lchoose <- function(N, k) {
    result <- vector(length=max(length(N), length(k)), mode="numeric")
    index <- which(N - k + 1 > 0)
    if (length(index) == 0) {
        result[] <- -Inf }
    else {
        result[index] <- (lgamma(N + 1) - lgamma(k + 1))[index] - lgamma((N - k + 1)[index])
        result[-index] <- -Inf
    }
    result
}



#    Copyright (C) 2016 University of Southern California and
#             Chao Deng and Andrew D. Smith and Timothy Daley
#
#    Authors: Chao Deng
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

## continued fraction approximant to a power series
## QD algorithm
## input: coefficients of the power series; begin with the constant
## mt: the maximum number of terms used in the power series
ps2cfa <- function(coef, mt) {
    index <- which(coef == 0)
    if (length(index) == 0) {
        mt <- min(mt, length(coef))
    } else {
        mt <- min(mt, index[1] - 1)
    }
    if (mt == 1) {
        return(coef[1])
    }
    qd.table <- matrix(data=0, nrow=mt, ncol=mt)
    ## initialize the table
    ## the first column is 0
    qd.table[1:(mt-1), 2] <- coef[2:mt] / coef[1:(mt-1)]
    if (mt == 2) {
        return(c(coef[1], -qd.table[1, 2]))
    }
    ## two types of columns e or q
    for (i in 3:mt) {
        n <- mt + 1 - i
        if (i %% 2 == 1) {
            ## number of entries in the column
            qd.table[1:n, i] <- qd.table[2:(n+1), i-1] - qd.table[1:n, i-1] + qd.table[2:(n+1), i-2]
            if (!is.finite(qd.table[1, i]) || qd.table[1, i] == 0) 
                return(c(coef[1], -qd.table[1, 2:(i-1)]))
        } else {
            qd.table[1:n, i] <- qd.table[2:(n+1), i-1] / qd.table[1:n, i-1] * qd.table[2:(n+1), i-2]
            if (!is.finite(qd.table[1, i]) || qd.table[1, i] == 0)
                return(c(coef[1], -qd.table[1, 2:(i-1)]))
        }
    }
    return(c(coef[1], -qd.table[1, 2:mt]))
}


## convert truncated continued fraction to a series of rational functions
## output two sets: A for numerators and B for denumerators
cfa2rf <- function(CF) {
    ## A, B are sets of polynomials based on recursive formula
    A <- list()
    B <- list()
    if (length(CF) < 2) {
        return(polynomial(CF))
    }
    A[[1]] <- polynomial(CF[[1]])
    A[[2]] <- polynomial(CF[[1]])
    B[[1]] <- polynomial(1)
    B[[2]] <- polynomial(c(1, CF[[2]]))
    if (length(CF) == 2) {
        return(list(A=A, B=B))
    }
    for (i in 3:length(CF)) {
        A[[i]] <- A[[i-1]] + polynomial(c(0, CF[[i]])) * A[[i-2]]
        B[[i]] <- B[[i-1]] + polynomial(c(0, CF[[i]])) * B[[i-2]]
    }
    return(list(A=A, B=B))
}

## Pad\'{e} approximant by picking out the numerator and the denominator
## input: two sets of polynomials for numerators and denominators
##        the degree m
## output: rational function approximant or Pad\'{e} approximant
rf2rfa <- function(RF, m) {
    return(polylist(RF$A[[m]], RF$B[[m]]))
}

## discriminant of the quadratic polynomial, which is
## the denominator of the discovery rate at m = 2
discriminant <- function(n) {
    if (max(n[, 1]) < 3) {
        return(NULL)
    }
    n[, 2] <- as.numeric(n[, 2])
    S1 <- sum(n[, 2])
    if (length(which(n[, 1] == 1))) {
        S2 <- S1 - n[which(n[, 1] == 1), 2]
    } else {
        S2 <- S1
    }
    if (length(which(n[, 1] == 2))) {
        S3 <- S2 - n[which(n[, 1] == 2), 2]
    } else {
        S3 <- S2
    }
    if (length(which(n[, 1] == 3))) {
        S4 <- S3 - n[which(n[, 1] == 3), 2]
    } else {
        S4 <- S3
    }
    a <- S2*S4 - S3^2
    b <- S1*S4 - S2*S3
    c <- S1*S3 - S2^2
    return((b / a)^2 - 4 * (c / a))
}