library(gtools)

rthetaB <- function(psi){
    
    return(theta)
}

rthetaB <- function(psi){
  K <- length(psi$gamma)
  pi <- as.vector(rdirichlet(1,psi$gamma))
  sig2 <- 1/rgamma(K,psi$lambda/2, psi$beta/2)
  mu <- rnorm(K,psi$alpha,sig2/psi$lambda)
  theta <- list(pi=pi,mu=mu,sig=sqrt(sig2))
  return(theta)
}

# bcp de variabilité avec le psi suivant :
psi <- list(gamma=c(1,1,1)*.1, alpha=c(1,1,1), lambda=c(1,1,1)*.5, beta=c(1,1,1))
rthetaB(psi)

# très peu de variabilité avec le psi suivant :
psi <- list(gamma=c(3,2,1)*100, alpha=c(1,1,1), lambda=c(1,1,1)*100, beta=c(1,1,1))
rthetaB(psi)

# pour des pi équilibrés : tous les gamma_k identiques (et une grande valeur pour gamma_k)
psi <- list(gamma=c(1,1,1)*30, alpha=c(1,1,1), lambda=c(1,1,1)*100, beta=c(1,1,1))
rthetaB(psi)

psi <- list(gamma=c(1,1,1)*30, alpha=c(1,1,1), lambda=c(1,1,1), beta=c(1,1,1)*100)
rthetaB(psi)

# 1: lambda grand et alpha distincts 
psi <- list(gamma=c(1,1,1)*100, alpha=c(1,10,20), lambda=c(1,1,1)*100, beta=c(1,1,1))
rthetaB(psi)
# 2: alpha identique, lambda petit et beta grand
psi <- list(gamma=c(1,1,1)*100, alpha=c(1,1,1), lambda=c(1,1,1)*100, beta=c(1,1,1)*100)
rthetaB(psi)

rnormmixB <- function(n,psi){

  return(list(obs=obs,theta=theta))
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

rnormmixB <- function(n,psi){
  theta <- rthetaB(psi)
  obs <- rnormmix(n,theta)
  return(list(obs=obs,theta=theta))
}

dnormmix <- function(x,theta){
  dens <- rep(0,length(x))
  K <- length(theta$pi)
  for (k in 1:K){
    dens <- dens + theta$pi[k]*dnorm(x,theta$mu[k],theta$sig[k])
  }                                         
  return(dens)
}

psi <- list(gamma=c(1,1,1)*3, alpha=c(0,5,15), lambda=c(1,1,1)*5, beta=c(1,1,1)*2)
res <- rnormmixB(100,psi)
res$theta
hist(res$obs,freq=F)
curve(dnormmix(x,res$theta),add=T,col='blue')

logaposteriori <- function(theta,x,psi){
    
    return(l)
}

logaposteriori <- function(theta,x,psi){
  vec <- log(dnormmix(x,theta))  
  vec[vec==-Inf] <- 0
  l <- sum(vec)
  l <- l + log(ddirichlet(theta$pi,psi$gamma))
  vec <- log(dnorm(theta$mu,psi$alpha,theta$sig/sqrt(psi$lambda)))
  vec[vec==-Inf] <- 0
  l <- l + sum(vec)
  vec <- log(dgamma(1/(theta$sig^2), (psi$lambda-1)/2,rate=psi$beta/2)/(theta$sig^4))
  vec[vec==-Inf] <- 0
  l <- l + sum(vec)
  return(l)
}

logaposteriori(res$theta,res$obs,psi)

n <- 200
psi <- list(gamma=c(3,3), alpha=c(0,4), lambda=c(4,4), beta=c(2,2))
set.seed(1711)
rthetaB(psi)
res <- rnormmixB(n,psi)
hist(res$obs)

mu.grille <- seq(-1.5,7,by=.1) 
G <- length(mu.grille)
L <- matrix(NA,nrow=G,ncol=G)
param <- res$theta
for (k in 1:G)
  for (l in 1:G){
    param$mu <- c(mu.grille[k],mu.grille[l]) 
    L[k,l] <- logaposteriori(param,res$obs,psi)     
  }  

image(mu.grille,mu.grille,L,col=grey(0:255/255))
contour(mu.grille,mu.grille,L,add=T,nlevels=20,col='red') 

rz <- function(theta,obs){

    return(z)
}

rz <- function(theta,obs){
  n <- length(obs)
  K <- length(theta$pi)
  P <- matrix(NA,n,K)
  for (k in 1:K)
    P[,k] <- theta$pi[k]*dnorm(obs,theta$mu[k],theta$sig[k])    
  z <- apply(P,1,function(row)if(sum(row)>0) sample(1:K,size=1,prob=row) else sample(1:K,size=1,prob=rep(1,K))) 
  return(z)
}

z <- rz(res$theta,res$obs)
z

rtheta <- function(z,obs,psi){
  
  return(theta)  
}

rtheta <- function(z,obs,psi){
  K <- length(psi$gamma)
  count_group <- sapply(1:K,function(k)sum(z==k))
  sum_group <- sapply(1:K,function(k)sum(obs[z==k]))
  sum2_group <- sapply(1:K,function(k)sum((obs[z==k])^2))
  
  pi <- as.vector(rdirichlet(1,psi$gamma+count_group))
  
  a <- (count_group+psi$lambda)/2
  b <- (sum2_group+psi$lambda*psi$alpha^2+psi$beta-(sum_group+psi$alpha*psi$lambda)^2/(count_group+psi$lambda))/2
  if (sum(b>0)==K){
    sig2 <- 1/rgamma(K,a,b)
  }else{ 
    sig2 <- 1/rgamma(K,(psi$lambda-1)/2, psi$beta/2)
  }
  m <- (sum_group+psi$lambda*psi$alpha)/(count_group+psi$lambda)
  s2 <- sig2/(count_group+psi$lambda)
  mu <- rnorm(K,m,sqrt(s2))
  
  theta <- list(pi=pi, mu=mu, sig=sqrt(sig2))
  return(theta)  
}

rtheta(z,res$obs,psi)

gibbs <- function(obs,theta.init,psi,R=1000,B=500){
 
  return(list(theta=theta,z=z))
}

gibbs <- function(obs,theta.init,psi,R=1000,B=500){
  theta <- vector("list",B+R+1)
  theta[[1]] <- theta.init
  z <- vector("list",B+R)
  for (r in 1:(R+B)){
    z[[r]] <- rz(theta[[r]],obs)
    theta[[r+1]] <- rtheta(z[[r]],obs,psi)
  }
  return(list(theta=theta[(B+2):(R+B+1)],z=z[(B+1):(R+B)]))
}

n <- 2000
psi <- list(gamma=c(3,3), alpha=c(0,4), lambda=c(4,4), beta=c(2,2))
set.seed(1711)
rthetaB(psi)
res <- rnormmixB(n,psi)

sim <- gibbs(res$obs,res$theta,psi)

mu <- sapply(sim$theta,function(theta) theta$mu  )
hist(mu[1,],freq=F,col='green',xlim=c(min(mu),max(mu)))
hist(mu[2,],freq=F,col=rgb(0,0,1,alpha=0.6),add=TRUE)

sig <- sapply(sim$theta,function(theta) theta$sig  )
hist(sig[1,],freq=F,col='green',xlim=c(min(sig),max(sig)))
hist(sig[2,],freq=F,col=rgb(0,0,1,alpha=0.6),add=TRUE)

pi <- sapply(sim$theta,function(theta) theta$pi)
hist(pi[1,],freq=F,col='green')

mu.grille <- seq(-1.5,7,by=.1) 
G <- length(mu.grille)
L <- matrix(NA,nrow=G,ncol=G)
param <- res$theta
for (k in 1:G)
  for (l in 1:G){
    param$mu <- c(mu.grille[k],mu.grille[l]) 
    L[k,l] <- logaposteriori(param,res$obs,psi)     
  }  

sim.mu <- sapply(sim$theta,function(elem) elem$mu) 

image(mu.grille,mu.grille,L,col=grey(0:255/255),xlab='mu_1',ylab='mu_2')
contour(mu.grille,mu.grille,L,add=T,nlevels=20,col='red') 
points(sim.mu[1,],sim.mu[2,],col='blue')

tracer.estim <- function(sim){
  
    return(estim)      
}

tracer.estim <- function(sim){
   sim.pi <- sapply(sim$theta,function(elem) elem$pi) 
   sim.mu <- sapply(sim$theta,function(elem) elem$mu) 
   sim.sig <- sapply(sim$theta,function(elem) elem$sig)     
   K <- nrow(sim.pi)
   R <- ncol(sim.pi)    
   par(mfrow=c(3,1))          
   plot(1:R,cumsum(sim.pi[1,])/(1:R),type='l',col=1,ylim=c(0,1),ylab='pi',xlab="nb d'itérations")    
   for(k in 2:K){
       lines(1:R,cumsum(sim.pi[k,])/(1:R),col=k)            
   }
   estim.mu <- apply(sim.mu,1,function(col)cumsum(col)/(1:R))    
   ymin <- min(estim.mu)
   ymax <- max(estim.mu)    
   plot(1:R,estim.mu[,1],type='l',col=1,ylim=c(ymin,ymax),ylab='mu',xlab="nb d'itérations")    
   for(k in 2:K){
       lines(1:R,estim.mu[,k],col=k)            
   }
   estim.sig <- apply(sim.sig,1,function(col)cumsum(col)/(1:R))    
   ymin <- min(estim.sig)
   ymax <- max(estim.sig)    
   plot(1:R,estim.sig[,1],type='l',col=1,ylim=c(ymin,ymax),ylab='sigma',xlab="nb d'itérations")    
   for(k in 2:K){
       lines(1:R,estim.sig[,k],col=k)            
   }
   estim=list(pi=rowMeans(sim.pi), mu=rowMeans(sim.mu), sig=rowMeans(sim.sig))
   return(estim)      
}

tracer.estim(sim)
res$theta

class.zi <- function(sim){
    
    return(class)
}

class.zi <- function(sim){
    K <- length(sim$theta[[1]]$pi)
    sim.z <- sapply(sim$z,function(elem)elem) # n x R
    class <- sapply(1:K,function(k) rowMeans(sim.z==k))
    return(class)
}

class <- class.zi(sim)
plot(res$obs,class[,1],xlab='x_i',ylab="probabilité d'être dans la classe 1")

theta.permute <- res$theta
theta.permute$pi <- theta.permute$pi[2:1]
theta.permute$mu <- theta.permute$mu[2:1]
theta.permute$sig <- theta.permute$sig[2:1]

sim <- gibbs(res$obs,theta.permute,psi)

sim.mu <- sapply(sim$theta,function(elem) elem$mu) 
image(mu.grille,mu.grille,L,col=grey(0:255/255),xlab='mu_1',ylab='mu_2')
contour(mu.grille,mu.grille,L,add=T,nlevels=20,col='red') 
points(sim.mu[1,],sim.mu[2,],col='blue')

theta <- res$theta
theta$mu <- c(3,2)  # entre les deux modes
sim <- gibbs(res$obs,theta,psi)
image(mu.grille,mu.grille,L,col=grey(0:255/255),xlab='mu_1',ylab='mu_2')
contour(mu.grille,mu.grille,L,add=T,nlevels=20,col='red') 
sim.mu <- sapply(sim$theta,function(elem) elem$mu) 
points(sim.mu[1,],sim.mu[2,],col='blue')   

class <- class.zi(sim)
plot(res$obs,class[,1],xlab='x_i',ylab="probabilité d'être dans la classe 1")

theta <- res$theta
theta$mu <- c(6,5) ##  logaposteriori très basse
sim <- gibbs(res$obs,theta,psi)
image(mu.grille,mu.grille,L,col=grey(0:255/255),xlab='mu_1',ylab='mu_2')
contour(mu.grille,mu.grille,L,add=T,nlevels=20,col='red') 
sim.mu <- sapply(sim$theta,function(elem) elem$mu) 
points(sim.mu[1,],sim.mu[2,],col='blue')

class <- class.zi(sim)
plot(res$obs,class[,1],xlab='x_i',ylab="probabilité d'être dans la classe 1")

load('donneesmelange.RData')

psi.ailes <- list(gamma=c(1/2,1/2), alpha= c(min(ailes),max(ailes)), lambda=c(2,2), beta=c(1,1)*var(ailes))
psi.ailes

theta.ailes <- list(pi=c(.5,.5),mu = c(min(ailes),max(ailes)),sig=sd(ailes)*c(1,1))
theta.ailes

gibbs.ailes <- gibbs(ailes,theta.ailes,psi.ailes,R=1000)

tracer.estim(gibbs.ailes)

class <- class.zi(gibbs.ailes)
plot(ailes,class[,1],xlab='x_i',ylab="probabilité d'être dans la classe 1")

psi.chlor <- list(gamma=c(1/2,1/2), alpha= c(min(chlor),max(chlor)), lambda=c(2,2), beta=c(1,1)*var(chlor))
psi.chlor

theta.chlor <- list(pi=c(.5,.5),mu = c(min(chlor),max(chlor)),sig=sd(chlor)*c(1,1))
theta.chlor

gibbs.chlor <- gibbs(chlor,theta.chlor,psi.chlor,R=5000)

tracer.estim(gibbs.chlor)

class <- class.zi(gibbs.chlor)
plot(chlor,class[,1],xlab='x_i',ylab="probabilité d'être dans la classe 1")
