library(boot)
boot(data, statistic, R, sim = "ordinary", stype = c("i", "f", "w"), 
     strata = rep(1,n), L = NULL, m = 0, weights = NULL, 
     ran.gen = function(d, p) d, mle = NULL, simple = FALSE, ...,
     parallel = c("no", "multicore", "snow"),
     ncpus = getOption("boot.ncpus", 1L), cl = NULL)
n <- 20
obs <- rnorm(n)
B <- 10
boot.mean <- function(obs, indice){
    moy.boot <- mean(obs[indice]) # permet à boot de choisir les données
    return(moy.boot)
}
repl.np <- boot(obs, boot.mean, B)
repl.np

names(repl.np)

repl.np$t

repl.np$t0

ran.gen.norm <- function(data, param){
  rnorm(length(data), param$mu, param$sig)
}
repl.param <- boot(obs, mean, B, sim='parametric', ran.gen=ran.gen.norm, 
                   mle=list(mu=mean(obs), sig=sd(obs)))
repl.param

dataexp <- rexp(20, 3)
T.exp.np <- boot(dataexp, function(obs, indice)1/mean(obs[indice]), R=1000)
T.exp.np

## génération d'un échnatillon bootstrap de loi exponentielle
ran.gen.exp <- function(data, param){
  rexp(length(data), param)
}



T.exp.param <- boot(dataexp, function(x)1/mean(x), 1000, sim='parametric',
                    ran.gen=ran.gen.exp, mle=1/mean(dataexp))
T.exp.param

boot.ci(boot.out, conf = 0.95, type = "all", 
        index = 1:min(2,length(boot.out$t0)), var.t0 = NULL, 
        var.t = NULL, t0 = NULL, t = NULL, L = NULL,
        h = function(t) t, hdot = function(t) rep(1,length(t)),
        hinv = function(t) t, ...)
B <- 1000
repl.np <- boot(obs, boot.mean, B)
icb.np <- boot.ci(repl.np)
icb.np

names(icb.np)

icb.np$normal

icb.np$basic

icb.np$percent

icb.np$bca

obs <- rnorm(20, 10, 2)       
boot.mean.var <- function(obs, indice){
    obs.boot <- obs[indice] 
    moy.boot <- mean(obs.boot) 
    var.boot <- var(replicate(50, mean(sample(obs.boot, replace=T))))
    return(c(moy.boot, var.boot))
}
repl.np.var <- boot(obs, boot.mean.var, R=1000)
repl.np.var

boot.ci(repl.np.var)$stud

## données exponentielles
dataexp <- rexp(20,3)
## répliques bootstrap non paramétrique
T.exp.np <- boot(dataexp, function(obs, indice)1/mean(obs[indice]), R=1000)
T.exp.np
## IC bootstrap non paramétrique    
confint.np <- boot.ci(T.exp.np, type=c("norm", "basic", "perc", "bca"))
confint.np

## IC studentisé, bootstrap non paramétrique
boot.estim.var.np <- function(obs, indice){
    obs.boot <- obs[indice]
    estim.boot <- 1/mean(obs.boot) 
    var.boot <- var(replicate(50, 1/mean(sample(obs.boot, replace=T))))
    return(c(estim.boot, var.boot))
}

# répliques bootstrap avec estimation de la variance
T.exp.var <- boot(dataexp, boot.estim.var.np, R=1000)
# ICB studentisé
confint.np <- boot.ci(T.exp.var)
confint.np

## ICB paramétriques
boot.estim.param <- function(obs){
    n <- length(obs)
    estim.boot <- 1/mean(obs) 
    var.boot <- var(replicate(50, 1/mean(rexp(n, estim.boot))))
    return(c(estim.boot, var.boot))
}

T.exp.param <- boot(dataexp, boot.estim.param, R=1000, sim='parametric',
                    ran.gen=ran.gen.exp, mle=1/mean(dataexp))
boot.ci(T.exp.param, type=c('norm', 'basic', 'perc'))

## ICB paramétrique type BCa 
boot.ci(T.exp.param, type=c('bca'))

## IC studentisé, bootstrap paramétrique
boot.estim.var.param <- function(obs){
    n <- length(obs)
    estim.boot <- 1/mean(obs) 
    var.boot <- var(replicate(50, 1/mean(rexp(n, estim.boot))))
    return(c(estim.boot, var.boot))
}

T.exp <- boot(dataexp, boot.estim.var.param, R=1000, sim='parametric',
              ran.gen=ran.gen.exp, mle=1/mean(dataexp))
T.exp
confint.param <- boot.ci(T.exp, type=c('norm', 'basic', 'stud', 'perc'))
confint.param

all.confint <- list(confint.param[4:7], confint.np[4:7])
all.confint

# theta = paramètre exponentiel
# n = taille d'échantillon
# M = nb de jeux de données simulés
# R = nb d'échantillons bootstrap
contains.theta <- function(vecIC, theta){
    L <- length(vecIC)
    is.theta.in <- (vecIC[L-1]<theta)&(vecIC[L]>theta)
    return(is.theta.in)
}

simul_icb_exp <- function(theta, n, M=100, R=100){
    prop <- matrix(0, 4, 2)
    colnames(prop) <- c('p','np')
    rownames(prop) <- c("norm", "basic", "stud", "perc")
    for (i in 1:M){
        obs <- rexp(n, theta)
        ## répliques bootstrap paramétriques avec estimations de variance
        T.param <- boot(obs, boot.estim.var.param, R, sim='parametric',
                        ran.gen=ran.gen.exp, mle=1/mean(obs))
        ## répliques bootstrap non paramétriques avec estimations de variance        
        T.np <- boot(obs, boot.estim.var.np, R)
        ## ICB paramétrique
        confint.param <- boot.ci(T.param, type=c('norm', 'basic', 'stud', 'perc'))
        ## ICB non paramétrique
        confint.np <- boot.ci(T.np,type=c("norm", "basic", "stud", "perc"))
        ## theta inclus dans ICB ...oui ou non 
        prop[,1] <- prop[,1] + sapply(confint.param[4:7], function(x) contains.theta(x, theta))
        prop[,2] <- prop[,2] + sapply(confint.np[4:7], function(x) contains.theta(x, theta))
    }
    return(prop/M)
}

simul_icb_exp(2, 10)

simul_icb_gamma <- function(gam, n, M=100, R=100){
    prop <- matrix(0, 4, 2)
    colnames(prop) <- c('p','np')
    rownames(prop) <- c("norm", "basic", "stud", "perc")
    for (i in 1:M){
        obs <- rgamma(n, gam$alpha, gam$theta)
        ## répliques bootstrap paramétriques avec estimations de variance
        T.param <- boot(obs, boot.estim.var.param, R, sim='parametric',
                        ran.gen=ran.gen.exp, mle=1/mean(obs))
        ## répliques bootstrap non paramétriques avec estimations de variance        
        T.np <- boot(obs, boot.estim.var.np, R)
        ## ICB paramétrique
        confint.param <- boot.ci(T.param, type=c('norm', 'basic', 'stud', 'perc'))
        ## ICB non paramétrique
        confint.np <- boot.ci(T.np, type=c("norm", "basic", "stud", "perc"))
        ## theta inclus dans ICB ?
        prop[,1] <- prop[,1] + sapply(confint.param[4:7], function(x) contains.theta(x, gam$theta/gam$alpha))
        prop[,2] <- prop[,2] + sapply(confint.np[4:7], function(x) contains.theta(x, gam$theta/gam$alpha))
    }
    return(prop/M)
}

simul_icb_gamma(list(alpha=.8, theta=2), 10)

simul_icb_gamma(list(alpha=.1, theta=2), 10)
