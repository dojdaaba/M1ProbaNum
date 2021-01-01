theta = list(pi=c(3,2,5)/10, mu=c(4,-5,0), sig=c(1.3,1,1.5))
theta

update.theta <- function(obs, theta){
  
  return(theta.new)
}

update.theta <- function(obs, theta){
  nb <- length(obs) 
  K  <- length(theta$mu) 
  
  # calcul des alpha_ji
  alpha <- matrix(NA,nrow=K, ncol=nb)
  for (j in 1:K){
    alpha[j,] <- theta$pi[j]*dnorm(obs,theta$mu[j],theta$sig[j])
  }
  p.theta <- apply(alpha,2,'sum')
  alpha   <- alpha/matrix(p.theta,byrow=TRUE,nrow=K,ncol=nb)
  
  # mise à jour des paramètres  
  pi.new  <- apply(alpha,1,'mean')
  mu.new  <- c(alpha%*%obs/nb/pi.new)
  mat     <- (matrix(rep(obs,each=K),nrow=K) -matrix(rep(mu.new,times=nb),nrow=K))^2
  sig.new <- sqrt(apply(alpha*mat,1,'mean')/pi.new)
  
  theta.new <-  list(pi =pi.new, mu=mu.new, sig=sig.new)
  return(theta.new)
}

rnormmix <- function(n,theta){
  K <- length(theta$pi)
  etiqu <- sample(1:K,size=n,replace=TRUE,prob=theta$pi)
  obs <- rep(NA,n)
  i <- 1
  for (k in 1:K){
    n.group <- sum(etiqu==k)
    if(n.group>0){
      obs[i:(i+n.group-1)] <- rnorm(n.group,theta$mu[k],theta$sig[k])
      i <- i + n.group
    }
  }        
  return(obs)
}

obs <- rnormmix(500,theta)
update.theta(obs,theta)

algoEM <- function(obs, theta.init){   

  resultat <- list(emv = theta.new, nb.iter = it)
  return(resultat)
}

dnormmix <- function(x,theta){
  dens <- rep(0,length(x))
  K <- length(theta$pi)
  for (k in 1:K){
    dens <- dens + theta$pi[k]*dnorm(x,theta$mu[k],theta$sig[k])
  }                                         
  return(dens)
}

lvraisnorm <- function(param,obs){    
  logvrais <- sum(log(dnormmix(obs,param)))
  return(logvrais)
}

algoEM <- function(obs, theta.init, R=200, epsilon=1e-3){   
  theta.old <- theta.init
  crit_arret <- FALSE
  log.vrais.old <- lvraisnorm(theta.init,obs)
  it <- 0
  while (!crit_arret && (it < R)){ 
    theta.new <- update.theta(obs, theta.old)
    log.vrais.new <- lvraisnorm(theta.new,obs)
    crit_arret <- (abs((log.vrais.new - log.vrais.old)/log.vrais.old) < epsilon)
    log.vrais.old <- log.vrais.new
    theta.old <- theta.new
    it <- it + 1
  }
  resultat <- list(emv = theta.new, nb.iter = it)
  return(resultat)
}

res <- algoEM(obs,theta)
res


ymax <- max(c(dnormmix(obs,res$emv),dnormmix(obs,theta)))
hist(obs,freq=FALSE,ylim=c(0,ymax))
curve(dnormmix(x,res$emv),add=TRUE,col='red')
curve(dnormmix(x,theta),add=TRUE,col='blue')

# permutation des composants et vraies valeurs
permut <- c(2,3,1)
theta1 <- list(pi=theta$pi[permut],mu=theta$mu[permut],sig=theta$sig[permut])
algoEM(obs,theta1)

# des mu_k très erronées
theta2 <- theta
theta2$mu <- c(10,20,30)
res2 <- algoEM(obs,theta2)
res2
hist(obs,freq=FALSE)
curve(dnormmix(x,res2$emv),add=TRUE,col='red')

# des mu_k très erronées
theta3 <- theta
theta3$mu <- c(-20,20,30)
res3 <- algoEM(obs,theta3)
res3
hist(obs,freq=FALSE)
curve(dnormmix(x,res3$emv),add=TRUE,col='red')

# des sig_k très erronées
theta4 <- theta
theta4$sig <-  c(10,10,10)  # c(.1,.1,.1), c(10,.1,10)
res4 <- algoEM(obs,theta4)
res4
hist(obs,freq=FALSE)
curve(dnormmix(x,res4$emv),add=TRUE,col='red')

# des pi_k erronés
theta5 <- theta
theta5$sig <-  c(1/3,1/3,1/3)  # c(.9, .1, .1), c(.4, .5, .1) c(.85, .01, .01)
res5 <- algoEM(obs,theta5)
res5
hist(obs,freq=FALSE)
curve(dnormmix(x,res5$emv),add=TRUE,col='red')

algoEM <- function(obs, K, theta.init=list(pi=rep(1/K,K),mu=sample(obs,size=K),sig=rep(sd(obs),K)), R=200, epsilon=1e-3){   
  # K - ordre du mélange
  theta.old <- theta.init
  crit_arret <- FALSE
  log.vrais.old <- lvraisnorm(theta.init,obs)
  it <- 0
  while (!crit_arret && (it < R)){ 
    theta.new <- update.theta(obs, theta.old)
    log.vrais.new <- lvraisnorm(theta.new,obs)
    crit_arret <- (abs((log.vrais.new - log.vrais.old)/log.vrais.old) < epsilon)
    log.vrais.old <- log.vrais.new
    theta.old <- theta.new
    it <- it + 1
  }
  resultat <- list(emv = theta.new, nb.iter = it)
  return(resultat)
}

res <- algoEM(obs,3)
res

theta

EMgauss <- function(obs,K,L){
    
    
    return(theta)
}

EMgauss <- function(obs,K,L=20){
  resEM <- vector("list",L)
  for (l in 1:L)
    resEM[[l]] <- algoEM(obs,K)$emv
  logvrais <- sapply(resEM, function(x) lvraisnorm(x,obs))
  best.run <- which.max(logvrais)    
  theta <- resEM[[best.run]]   
  return(theta)
}

EMgauss(obs,3)

load('donneesmelange.RData')
table(ailes)
table(chlor)

mu.ailes <- mean(ailes)
sig.ailes <- sd(ailes)
theta.ailes <- EMgauss(ailes,2)
theta.ailes
ymax <- max(c(dnorm(ailes,mu.ailes,sig.ailes),dnormmix(ailes,theta.ailes)))
hist(ailes,freq=FALSE,ylim=c(0,ymax))
curve(dnorm(x,mu.ailes,sig.ailes),add=TRUE,col='green')
curve(dnormmix(x,theta.ailes),add=TRUE,col='red')

set.seed(1984)
mu.chlor <- mean(chlor)
sig.chlor <- sd(chlor)
theta.chlor <- EMgauss(chlor,2,100)
theta.chlor
ymax <- max(c(dnorm(chlor,mu.chlor,sig.chlor),dnormmix(chlor,theta.chlor)))
hist(chlor,freq=FALSE,ylim=c(0,ymax))
curve(dnorm(x,mu.chlor,sig.chlor),add=TRUE,col='green')
curve(dnormmix(x,theta.chlor),add=TRUE,col='red')
