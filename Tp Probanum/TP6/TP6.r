plot_ly(x=x.grid, y=y.grid, z=z.values, type="surface")ggplot(df, aes(x=x, y=y, z=f_xy)) + geom_contour(bins = 60, aes(colour = ..level..))
log_vrais <- function(a,b,obs){
  n <- length(obs)
  log_vrais <- n*a*log(b)-n*log(gamma(a))+(a-1)*sum(log(obs))-b*sum(obs)
  return(log_vrais)
}

library(ggplot2)
library(plyr)

log_vrais3D <- function(obs, a.grid=seq(.1,5,by=.01), b.grid=seq(.1,5,by=.01)){
  L <- sapply(1:length(b.grid), function(i) return(log_vrais(a.grid,b.grid[i],obs)))
  A <- length(a.grid)
  B <- length(b.grid)
  df <- data.frame(a=rep(a.grid,times=B), b=rep(b.grid,each=A), log_vrais=c(L))  # création d'un df pour ggplot
  p <- ggplot(df, aes(a, b, z=log_vrais)) + geom_contour(bins = 60, aes(colour = ..level..))
  return(p)
}

obs <- rgamma(100,2,.8)
log_vrais3D(obs) 
update_newton <- function(param, obs){
    # param - liste de paramètres: param$a, param$b
    # obs - vecteur de données
    
    # calcul des nouvelles valeurs des paramètres :
    ...
    a.new <- ...
    b.new <- ...
    
    return(list(a=a.new,b=b.new))
}
update_newton <- function(param, obs){
  # param : liste des paramètres a et b
  # obs - vecteur de données
  n <- length(obs)
  gradient <- c(n*log(param$b)-n*digamma(param$a)+sum(log(obs)), n*param$a/param$b-sum(obs))
  w <- trigamma(param$a)
  detH <- 1/n/(1-param$a*w)
  # calcul des nouvelles valeurs des paramètres :
  a.new <- param$a - detH*sum(c(param$a,param$b)*gradient)
  b.new <- param$b - detH*sum(c(param$b,param$b^2*w)*gradient)
  return(list(a=a.new,b=b.new))
}
emv_gamma <- function(obs, param){
    # obs - vecteur de données
    # param - liste de paramètres: param$a, param$b
    
    param.old <- param
    not.converged <- TRUE
    ...
    while(not.converged){
        # mise à jour des paramètres :
        param.new <- update_newton(param,obs)
        # vérifier le critère d'arrêt :
        not.converged <- ...
        param.old <- param.new
        ...
    }
    
    return(list(param=param, nb.iter=nb.iter))    
}
emv_gamma <- function(obs,param){
  # obs - vecteur de données
  # param - liste de paramètres: param$a, param$b
  epsilon <- 10^{-3}
  param.old <- param
  log_vrais.old <- log_vrais(param.old$a,param.old$b,obs)
  not.converged <- TRUE
  nb.iter <- 0
  while (not.converged){
    # mise à jour des paramètres :
    param.new <- update_newton(param.old,obs)
    # vérifier le critère d'arrêt :
    log_vrais.new <- log_vrais(param.new$a,param.new$b,obs)
    not.converged <- (abs((log_vrais.new-log_vrais.old)/log_vrais.new)>epsilon)
    param.old <- param.new   
    log_vrais.old <- log_vrais.new
    nb.iter <- nb.iter + 1
  }
  return(list(param=param.new,nb.iter=nb.iter))  
}

emv_gamma <- function(obs, param){
  # obs - vecteur de données
  # param - liste de paramètres: param$a, param$b
  epsilon <- 10^{-3}
  nb.iter.max <- 100  
  param.old <- param
  log_vrais.old <- log_vrais(param.old$a, param.old$b,obs)
  not.converged <- TRUE
  nb.iter <- 0
  while ((not.converged) & (nb.iter <= nb.iter.max)){  #modification
    # mise à jour des paramètres :
    param.new <- update_newton(param.old,obs)
    # vérifier le critère d'arrêt :
    log_vrais.new <- log_vrais(param.new$a,param.new$b,obs)
    not.converged <- (abs((log_vrais.new-log_vrais.old)/log_vrais.new)>epsilon)
    param.old <- param.new
    log_vrais.old <- log_vrais.new
    nb.iter <- nb.iter + 1
  }
  return(list(param=param.new, nb.iter=nb.iter))  
}

emv_gamma(rgamma(100,3,3), list(a=3,b=3))

emv_gamma(rgamma(200,0.2,4), list(a=1,b=1))

emm_gamma <- function(obs){
  n <- length(obs)
  b <- mean(obs)/(var(obs)*(n-1)/n)
  a <- b*mean(obs)
  return(list(a=a, b=b))
}

emv_gamma <- function(obs, param=emm_gamma(obs)){
  # obs - vecteur de données
  # param - liste de paramètres: param$a, param$b
  epsilon <- 10^{-3}
  nb.iter.max <- 100  
  param.old <- param
  log_vrais.old <- log_vrais(param.old$a, param.old$b,obs)
  not.converged <- TRUE
  nb.iter <- 0
  while ((not.converged) & (nb.iter <= nb.iter.max)){
    # mise à jour des paramètres :
    param.new <- update_newton(param.old, obs)
    if (sum(c(param.new)>0)==2){
      # vérifier le critère d'arrêt :
      log_vrais.new <- log_vrais(param.new$a, param.new$b,obs)
      not.converged <- (abs((log_vrais.new-log_vrais.old)/log_vrais.new)>epsilon)
      param.old <- param.new
      log_vrais.old <- log_vrais.new
    }
    else{
      not.converged <- FALSE
      }
    nb.iter <- nb.iter + 1  
  }
  return(list(param=param.new, nb.iter=nb.iter))  
}

# voilà le jeu de données de tout à l'heure qui ne marchait pas:
set.seed(12)
obs <- rgamma(200,.2,4)
emv_gamma(obs)

emv_gamma(rgamma(50,2,3))

library(datasets)
param <- emv_gamma(lynx)$param
param
hist(lynx, breaks=30, freq=FALSE)
curve(dgamma(x, param$a, param$b), add=TRUE,col='red')

comp_risk <- function(a, b, n, R=1000, seed=17){
    
  # a,b - paramètres de la loi Gamma
  # n - taille d'échantillon
  # R - nb d'échantillon à simuler
  
    
    # evaluation des 2 estimateurs 
  set.seed(seed) # pour utiliser les mêmes données pour les deux estimateurs
  emv <- replicate(R,as.numeric(emv_gamma(rgamma(n, a, b))$param))  
  risk_emv <- sum((emv - matrix(c(a, b), nrow=2, ncol=R))^2)/R
  
  set.seed(seed) 
  emm <- replicate(R,as.numeric(emm_gamma(rgamma(n, a, b))))  
  risk_emm <- sum((emm - matrix(c(a, b), nrow=2, ncol=R))^2)/R
  
  return(list(risk_emm=risk_emm, risk_emv=risk_emv))
}

comp_risk(1, 2, 200)

comp_risk(1, 2, 2000)
