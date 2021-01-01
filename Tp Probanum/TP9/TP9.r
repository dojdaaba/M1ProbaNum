theta = list(pi=c(.2,.2,.6), mu=c(2,-5,2), sig=c(2,1,.1))
theta

dnormmix <- function(x, theta){
    
    return(dens)
}

dnormmix <- function(x, theta){
    dens <- rep(0, length(x))
    K <- length(theta$pi)
    for (k in 1:K){
        dens <- dens + theta$pi[k]*dnorm(x, theta$mu[k], theta$sig[k])
    }                                         
    return(dens)
}

par(mfrow=c(2,3))
# bimodal symétrique, proche
theta <- list(pi=c(.5,.5), mu = c(-1,1), sig= c(2/3,2/3))
curve(dnormmix(x,theta),-3,3,lwd=2,ylab='',xlab='')

# trimodal, séparé
theta <- list(pi=c(.4,.3,.3), mu = c(-1,3,7), sig= c(1.2,2/3,1.4))
curve(dnormmix(x,theta),-3,10,lwd=2,ylab='',xlab='')

# pique, coin
theta <- list(pi=c(2/3,1/3), mu = c(0,0), sig= c(1,1/10))
curve(dnormmix(x,theta),-3,3,lwd=2,ylab='',xlab='')

# bimodal, proche, asymétrique
theta <- list(pi=c(3/4,1/4), mu = c(0,3/2), sig= c(1,1/3))
curve(dnormmix(x,theta),-3,3,lwd=2,ylab='',xlab='')

# asymétrique
theta <- list(pi=c(.4,.15,.45), mu = c(0,1/2,13/15), sig= c(1,2/3,5/9))
curve(dnormmix(x,theta),-3,3,lwd=2,ylab='',xlab='')

# peigne
theta <- list( pi=c(1/2,rep(1/10,5)), mu = c(0,(0:4)/2-1), sig= c(1,rep(1/10,5)) )
x <- seq(-2.5,2.5,by=.01)
plot(x,dnormmix(x,theta),lwd=2,ylab='',xlab='',type='l')

dunifmix <- function(x, theta){
    # theta : $pi - poids de composant
    #         $lambda - paramètres uniformes
    dens <- rep(0, length(x))
    K <- length(theta$pi)
    for (k in 1:K){
        dens <- dens + theta$pi[k]*dunif(x, 0, theta$lambda[k])
    }                                         
    return(dens)
}

par(mfrow=c(1,3))
theta1 = list(pi=c(1,1)/2, lambda=c(2,4))
curve(dunifmix(x,theta1),-.5,5)
theta1 = list(pi=c(1,2,3)/6, lambda=c(1,2,3))
curve(dunifmix(x,theta1),-.5,5)
theta2 = list(pi=rep(1/50,50), lambda=(1:50)/50)
curve(dunifmix(x,theta2),-.5,1.5)

dgammamix <- function(x, theta){
    # theta : $pi - poids de composant
    #         $alpha, $beta - paramètres de la loi Gamma
    dens <- rep(0, length(x))
    K <- length(theta$pi)
    for (k in 1:K){
        dens <- dens + theta$pi[k]*dgamma(x, theta$alpha[k], theta$beta[k])
    }                                         
    return(dens)
}

par(mfrow=c(2,2))
theta1 = list(pi=c(1,1)/2,alpha=c(2,4),beta=c(2,1))
curve(dgammamix(x, theta1),-.5,15)

theta2 = list(pi=c(1,4)/5,alpha=c(2,20),beta=c(2,2))
curve(dgammamix(x, theta2),-.5,20)

theta3 = list(pi=rep(1/5,5),alpha= c(1,9,22,38,60) ,beta=rep(2,5))
curve(dgammamix(x, theta3),-.5,45)

theta4 = list(pi=c(1,4,6,8,11)/30,alpha= c(1,9,22,38,60) ,beta=rep(2,5))
curve(dgammamix(x, theta4),-.5,45)


rnormmix <- function(n, theta){

    return(obs)
}

rnormmix <- function(n, theta){
  K <- length(theta$pi)
  etiqu <- sample(1:K, size=n, replace=TRUE, prob=theta$pi)
  obs <- rep(NA,n)
  i <- 1
  for (k in 1:K){
    n.group <- sum(etiqu==k)
    if(n.group>0){
      obs[i:(i+n.group-1)] <- rnorm(n.group, theta$mu[k], theta$sig[k])
      i <- i + n.group
    }
  }        
  return(obs)
}

theta <- list(pi=c(.5,.5), mu = c(-1,1), sig= c(2/3,2/3))
obs <- rnormmix(1000, theta)
hist(obs, freq=F)
curve(dnormmix(x,theta), add=T)

lvraisnorm <- function(param, obs){
    
    return(logvrais)
}

lvraisnorm <- function(param, obs){    
    logvrais <- sum(log(dnormmix(obs, param)))
    return(logvrais)
}

lvraisnorm(theta, obs)
image(x, y, z, zlim, xlim, ylim, col = heat.colors(12),
      add = FALSE, xaxs = "i", yaxs = "i", xlab, ylab,
      breaks, oldstyle = FALSE, useRaster, ...)
x <- seq(0, 2, by=.03)
y <- seq(0, 1, by=.03)
M <- matrix(NA, nrow=length(x), ncol=length(y))
for (k in 1:length(x))
  M[k,] <- y^x[k]
image(x, y, M, col=grey(0:255/255))
contour(x,y,M,add=T,nlevels=10,col='yellow')

#contour(x, y, M, add=T, nlevels=10, col='yellow')

theta <- list(pi=c(1,1)/2, mu=c(0,2.5), sig=c(1,1))
obs <- rnormmix(1000, theta)
mu1 <- seq(-1, 4, by=.02)
mu2 <- seq(-2, 4, by=.02)
L <- matrix(NA, nrow=length(mu1), ncol=length(mu2))
param <- theta
for (k in 1:length(mu1))
    for (l in 1:length(mu2)){
        param$mu <- c(mu1[k], mu2[l]) 
        L[k,l] <- lvraisnorm(param, obs)
    }


image(mu1, mu2, L, col=grey(0:255/255))
contour(mu1, mu2, L, add=T, nlevels=50, col='red')

theta <- list(pi=c(3,7)/10, mu=c(0,2.5), sig=c(1,1))
obs <- rnormmix(1000, theta)
mu1 <- seq(-1, 4, by=.02)
mu2 <- seq(-2, 4, by=.02)
L <- matrix(NA, nrow=length(mu1), ncol=length(mu2))
param <- theta
for (k in 1:length(mu1))
    for (l in 1:length(mu2)){
        param$mu <- c(mu1[k], mu2[l]) 
        L[k,l] <- lvraisnorm(param, obs)
    }

image(mu1, mu2, L, col=grey(0:255/255))
contour(mu1, mu2, L, add=T, nlevels=50, col='red')
