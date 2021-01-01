nsim <- 500
x <- runif(nsim, min = -4, max = 8)
plot(x, ylim = c(-4,8), ylab = NA, pch = 3, cex=.5, col="grey", main = "Illustration TCL, loi uniforme")

abline(h = 2, col = "blue")                     # la ligne bleue 
m <- cumsum(x) / seq_along(x)                   # le vecteur m contient les moyennes cumulées 
                                                # on peut aussi diviser par (1:numsim)
lines(m, col = "red")                           

v <- sapply(2:nsim, function(k)               # on calcule la variance pour $n> 2$
            sum((x[1:k] - m[k])^2) / (k-1)) 
v <- c(NA, v)            # NA pour n = 1
head(v)
# autre approche équivalente
v <- sapply(1:nsim, function(k) var(x[1:k]))  # on calcule directement la variance avec la fonction var
head(v)
e <- 1.96 * sqrt(v) / sqrt(1:nsim)            # taille de l'IC à 95%
lines(m+e)
lines(m-e)

numsim <- 5000
x <- rcauchy(numsim)
plot(x, ylim = c(-10,10), ylab = NA, pch = 3, cex=.5, col="grey", main = "Illustration TCL, loi de Cauchy")

abline(h = 0, col = "blue")                    
m <- cumsum(x) / seq_along(x)                  
lines(m, col = "red")                          
v <- sapply(1:numsim, function(k) var(x[1:k]))   
e <- 1.96 * sqrt(v) / sqrt(1:numsim)          
lines(m+e)
lines(m-e)

n <- 1e6
monte_carlo <- function(beta) {
    x <- rnorm(n)
    y <- exp(beta * x)
    m <- mean(y)
    v <- var(y)
    e <- 1.96*sqrt(var(y)/length(y))
    return(c(m, m-e, m+e, v))
}

beta <- c(0.2, 0.5, 1, 2, 3, 5)
r <- sapply(beta, monte_carlo)
colnames(r) <- beta
rownames(r) <- c("Mean", "Lower", "Upper", "Var")
                           
rbind(r)     

d <- 10
T <- 1
r <- 0.04
S_0 <- rep(100, d)
sigma <- 1:d/(2*d)
rho <- 0.5
Sigma <- matrix(0.2, 10, 10) + diag(1-0.2, 10, 10)
L <- t(chol(Sigma))

phi <- function(G, K) {
    S_T <- S_0 * exp((r - 0.5 * sigma^2) * T + sigma * sqrt(T) * (L %*% G))
    return(max(sum(S_T)/d-K, 0))   
}

monte_carlo <- function(n, X, Ks) {
    result <- matrix(nrow = 0, ncol = 5)
    colnames(result) <- c("K", "Mean", "Var", "Lower", "Upper")
    for (K in Ks) {
        x <- replicate(n, X(K))
        mv <- c(mean(x), var(x))
        e_n <- 1.96 * sqrt(mv[2]) / sqrt(n) 
        result <- rbind(result, c(K, mv, mv[1]-e_n, mv[1]+e_n))  
    }
    return(result)
}

n <- 1e5
Ks <- c(80, 90, 100, 110, 120)
X <- function(K) { return(phi(rnorm(d, 0, 1), K)) }
r_mc <- monte_carlo(n, X, Ks)
r_mc

n <- 1e5
Ks <- c(80, 90, 100, 110, 120)
X_antith <- function(K) {
    G <- rnorm(d, 0, 1) 
    return(0.5 * (phi(G, K) + phi(-G, K)))
}
r_antith <- monte_carlo(n, X_antith, Ks)
r_antith

ratios <- cbind(r_antith[,"K"], r_mc[,"Var"] / r_antith[,"Var"])              # le ratio des variances
colnames(ratios) <- c("K", "Var MC / Var Antith")
ratios

bs_call <- function(m, s2, K) {
    d1 <- (m - log(K)) / sqrt(s2)
    return(exp(m + s2/2) * pnorm(d1 + sqrt(s2)) - K * pnorm(d1))
}

esp_gZ <- function(K) {
    meanZ <- (T/d) * sum(r - (sigma^2) / 2)
    varZ <- 0
    for (i in 2:d) {
        for (j in 1:(i-1)) {
            varZ <- varZ + 2 * rho * sigma[i]*sigma[j]
        }
    }
    varZ <- T / (d^2) * (sum(sigma^2) + varZ)
    return(S_0[1] * bs_call(meanZ, varZ, K / S_0[1]))
}

esp_gZ(Ks)

psi <- function(G, K) {
    Z <- sum((r-0.5*sigma^2)*T + sigma*sqrt(T)*L %*% G) / d
    return(max(S_0 * exp(Z) - K, 0))
}

n <- 1e5
Ks <- c(80, 90, 100, 110, 120)
X_vc <- function(K) {
    G <- rnorm(d, 0, 1) 
    return(phi(G, K) - psi(G, K))
}
r_vc <- monte_carlo(n, X_vc, Ks)

corrections <- esp_gZ(Ks)
cols <- c("Mean", "Lower", "Upper")
r_vc[,cols] <- r_vc[,cols] + corrections
r_vc

ratios <- cbind(r_vc[,"K"], r_mc[,"Var"] / r_vc[,"Var"])              
colnames(ratios) <- c("K", "Var MC / Var VC")
ratios



v <- rep(1/sqrt(d), d)
A <- diag(d) - v %*% t(v)
Nstrates <- 10
pi <- 1 / Nstrates
Z <- function(i) {
    G <- rnorm(d, 0, 1)
    U <- runif(1, 0, 1)
    Ui <- qnorm((i-1)/Nstrates + U / Nstrates, 0, 1)
    return(Ui * v + A %*% G)
}

n <- 10e4
ni <- n * pi
mcs <- numeric(Nstrates)
vars <- numeric(Nstrates)
K <- 100
for (i in 1:Nstrates) {
    xi <- replicate(ni, phi(Z(i), K))
    mcs[i] <- mean(xi) * pi
    vars[i] <- var(xi) * pi
}
mcs
sum(mcs)
sum(vars)
