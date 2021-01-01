aircond <- c(487, 18, 85, 91, 5, 130, 230, 43, 98, 7, 100, 3)

hist(aircond, freq=FALSE, breaks=8)
curve(dexp(x, 1/mean(aircond)), add=TRUE)

obs1.boot <- sample(aircond, replace=TRUE)
obs2.boot <- sample(aircond, replace=TRUE)
obs3.boot <- sample(aircond, replace=TRUE)
plot(ecdf(obs1.boot), xlim=c(0, max(aircond)), col=2, main='Echantillons bootstrap non paramétrique')
plot(ecdf(obs2.boot), add=TRUE, col=3)
plot(ecdf(obs3.boot), add=TRUE, col=4)
plot(ecdf(aircond), add=TRUE, col=1)
legend(300, .6, lty=1, legend='aircond')

n <- length(aircond)
theta.hat <- 1/mean(aircond)
obs1.boot <- rexp(n, theta.hat)
obs2.boot <- rexp(n, theta.hat)
obs3.boot <- rexp(n, theta.hat)
plot(ecdf(obs1.boot), xlim=c(0, max(aircond)), col=2, main='Echantillons bootstrap paramétrique')
plot(ecdf(obs2.boot), add=TRUE, col=3)
plot(ecdf(obs3.boot), add=TRUE, col=4)
plot(ecdf(aircond), add=TRUE, col=1)
curve(pexp(x, theta.hat), add=TRUE,col=5)
legend(300, .6, lty=1, legend='aircond')

boot.np <- function(obs, R=1000){
  return(replicate(R, 1/mean(sample(obs, replace=TRUE))))
}

boot.exp <- function(obs, R=1000){
  n <- length(obs)
  theta.hat <- 1/mean(obs)  
  return(replicate(R, 1/mean(rexp(n, theta.hat))))
}

estim.T <- function(theta, n, R=1000){
    return(replicate(R, 1/mean(rexp(n, theta))))
}

estim.np <- boot.np(aircond, 1000)
estim.exp <- boot.exp(aircond, 1000)

par(mfrow=c(1, 2))
hist(estim.np, freq=FALSE, main='Bootstrap nonparamétrique', xlab='répliques bootstrap')
hist(estim.exp, freq=FALSE, main='Bootstrap paramétrique', xlab='répliques bootstrap')

centreduit <- function(vec)
    return((vec-mean(vec))/sd(vec))

par(mfrow=c(1,2))
qqnorm(centreduit(estim.np), ylab='répliques non paramétriques')
abline(0,1)
qqnorm(centreduit(estim.exp), ylab='répliques paramétriques')
abline(0,1)

plot(ecdf(estim.np), col='green', main='Répliques bootstrap')
plot(ecdf(estim.exp), col='blue', add=TRUE)
legend(.015, .5, legend=c('nonparamétrique', 'paramétrique'), col=c('green', 'blue'), lty=1)

qqplot(estim.np, estim.exp, xlab='répliques non paramétriques', ylab='répliques paramétriques', main='QQ-plot')
abline(0, 1)

qqplot(centreduit(estim.np), centreduit(estim.exp), xlab='répliques non paramétriques', ylab='répliques paramétriques', main='QQ-plot')
abline(0, 1)

n <- length(aircond)
theta <- 1/mean(aircond)
T <- estim.T(theta, n, 100000)
par(mfrow=c(1, 2))
qqplot(T, estim.np, main='Bootstrap non paramétrique')
abline(0,1)
qqplot(T, estim.exp, main='Bootstrap paramétrique')
abline(0,1)

simul_min <- function(theta, n, R=1000){
  estim <- replicate(R, min(runif(n, theta, theta+1)))
  Ln <- n*(estim-theta)  
  par(mfrow=c(1, 2))
  hist(Ln, freq=FALSE)
  curve(dexp(x, 1),col=2, add=TRUE)
  qqplot(qexp((1:n)/n, 1), sort(Ln), xlab='Loi exponentielle')
  abline(0, 1)
  return(Ln)
}

Ln <- simul_min(2, 50, 1000)

boot_min <- function(obs, R=1000){
  estim <- replicate(R, min(sample(obs, replace=TRUE)))
  n <- length(obs)
  theta <- min(obs)
  Ln <- n*(estim-theta)  
  par(mfrow=c(1, 2))
  hist(Ln, freq=FALSE)
  curve(dexp(x, 1), col=2, add=TRUE)
  qqplot(qexp((1:n)/n, 1), Ln, xlab='Loi exponentielle')
  abline(0, 1)
  return(Ln)
}

ln <- boot_min(runif(10000, 2, 3), 1000)

ln <- boot_min(runif(100, 2, 3), 1000)
unique(ln)

bootparam_min <- function(obs, R=1000){
  n <- length(obs)
  theta <- min(obs)
  estim <- replicate(R, min(runif(n, theta, theta+1)) )
  Ln <- n*(estim-theta)  
  par(mfrow=c(1, 2))
  hist(Ln, freq=FALSE)
  curve(dexp(x, 1), col=2, add=TRUE)
  qqplot(qexp((1:n)/n, 1), Ln, xlab='Loi exponentielle')
  abline(0, 1)
  return(Ln)
}

ln <- bootparam_min(runif(50, 2, 3), 1000)

biais.var.boot <- function(obs, B, param=FALSE){
    if (param)
        estim.boot <- boot.exp(obs, B)
    else
        estim.boot <- boot.np(obs, B)
    b <- mean(estim.boot)- 1/mean(obs)
    v <- var(estim.boot)    
    return(list(biais = b, var = v))
}

obs <- rexp(15, .1)

biais.var.boot(obs, 300, param=TRUE)

biais.var.boot(obs, 300, param=TRUE)

icb.norm <- function(obs, alpha=.05, B=300, param=FALSE){
    bv <- biais.var.boot(obs, B, param)
    theta <- 1/mean(obs)
    return(theta-bv$biais -sqrt(bv$var)*qnorm(1-alpha/2)*c(1, -1))
}

icb.norm(obs, param=FALSE)

niveau.icbnorm <- function(theta, n, M, alpha=0.05, param=FALSE){
    prop <- 0
    for (i in 1:M){
        obs <- rexp(n, theta)
        icb <- icb.norm(obs, alpha, param=param)
        prop <- prop + ((icb[1]<theta)&(theta<icb[2]))
    }
    return(prop/M)
}

niveau.icbnorm(2, 10, 1000, param=TRUE)

niveau.icbnorm(1/mean(aircond), 12, 1000, param=TRUE)

niveau.icbnorm(.02, 10, 1000, param=FALSE)

icb.basic <- function(obs, alpha=.05, R=1000, param=FALSE){
    gamma <- c(1-alpha/2, alpha/2)
    theta <- 1/mean(obs)
    if (param)
        estim.boot <- boot.exp(obs, R)
    else
        estim.boot <- boot.np(obs, R)
   quant <- as.numeric(quantile(estim.boot, gamma))
    return(2*theta-quant)
}

icb.basic(obs, param=TRUE)

niveau.icbbasic <- function(theta, n, M, alpha=0.05, param=FALSE){
    prop <- 0
    for (i in 1:M){
        obs <- rexp(n, theta)
        icb <- icb.basic(obs, alpha, param=param)
        prop <- prop + ((icb[1]<theta)&(theta<icb[2]))
    }
    return(prop/M)
}

niveau.icbbasic(.2, 10, 1000, param=TRUE)

niveau.icbbasic(.2, 50, 1000, param=FALSE)

niveau.icbbasic(.2, 50, 1000, param=FALSE)

niveau.icbbasic(.2, 50, 1000, param=TRUE)

icb.stud <- function(obs, alpha=.05, R=300, B=50, param=FALSE){
    gamma <- c(1-alpha/2, alpha/2)
    theta <- 1/mean(obs)
    n <- length(obs)
    estim.boot <- rep(NA, R)
    Z <- rep(NA, R)
    for (i in 1:R){
        if (param)
            obs.boot <- rexp(n, theta)
        else    
            obs.boot <- sample(obs, replace=TRUE)            
        estim.boot[i] <- 1/mean(obs.boot)
        var.boot <- biais.var.boot(obs.boot,B)$var
        Z[i] <- (estim.boot[i]-theta)/sqrt(var.boot)
    }
    quant <- as.numeric(quantile(Z, gamma))
    var.boot <- var(estim.boot)
    return(theta-sqrt(var.boot)*quant)
}

icb.stud(obs, param=TRUE, R=1000)

niveau.icbstud <- function(theta, n, M, alpha=0.05, param=FALSE){
    prop <- 0
    for (i in 1:M){
        obs <- rexp(n, theta)
        icb <- icb.stud(obs, alpha, param=param)
        prop <- prop + ((icb[1]<theta)&(theta<icb[2]))
    }
    return(prop/M)
}

niveau.icbstud(2, 10, 100, param=TRUE)

niveau.icbstud(2, 10, 100, param=FALSE)

niveau.icbstud(.2, 30, 100, param=TRUE)

niveau.icbstud(.2, 30, 100, param=FALSE)


icb.perc <- function(obs, alpha=.05, R=1000, param=FALSE){
    gamma <- c(alpha/2, 1-alpha/2)
    if (param)
        estim.boot <- boot.exp(obs, R)
    else
        estim.boot <- boot.np(obs, R)
    quant <- as.numeric(quantile(estim.boot, gamma))
    return(quant)
}

icb.perc(obs)

niveau.icbperc <- function(theta, n, M, alpha=0.05, param=FALSE){
    prop <- 0
    for (i in 1:M){
        obs <- rexp(n, theta)
        icb <- icb.perc(obs, alpha, param=param)
        prop <- prop + ((icb[1]<theta)&(theta<icb[2]))
    }
    return(prop/M)
}

niveau.icbperc(2, 10, 1000, param=FALSE)

niveau.icbperc(2, 10, 1000, param=TRUE)

niveau.icbperc(2, 30, 1000, param=FALSE)

niveau.icbperc(2, 30, 1000, param=TRUE)

icb.norm(aircond)

icb.basic(aircond)

icb.stud(aircond)

icb.perc(aircond)
