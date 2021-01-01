lambda <- 10
mu <- 0.1
sigma <- 0.3
S <- function() {
  N <- rpois(1, lambda)
  X <- rlnorm(N, mu, sigma)
  return(sum(X))
}
n <- 1e5
s <- replicate(n, S())
hist(s, breaks = 50, prob = T, main = paste("Histogramme d'un échantillon de S de taille n=", n))
K <- 20
abline(v = K, col = "red", lty = 1)
legend(22, 0.08, paste("seuil K = ", K), col = "red", lty = 1)

erreur_relative <- function(x, K) {
    p_n <- mean(x > K)
    v_n <- p_n - p_n^2
    e_n <- 1.96 * sqrt(v_n) / sqrt(length(x))
    return(e_n / p_n)
}
s <- replicate(1e5, S())
Ks <- seq(20,30,1)
er_K <- sapply(Ks, function(K) erreur_relative(s, K))
plot(Ks, er_K, t = "l", main="Erreur relative en fonction du seuil K")

m <- 1e4
monte_carlo <- function(K) {
    s <- replicate(m, S())
    e_n <- erreur_relative(s, K)
    print(e_n)
    while (e_n > 0.05) {
        new_sample <- replicate(m, S())
        s <- c(s, new_sample)
        e_n <- erreur_relative(s, K)
    }
    return(c(mean(s > K), e_n, length(s)))
}

Ks <- seq(20, 25, by = 1)
MC_K <- sapply(Ks, monte_carlo)
colnames(MC_K) <- Ks
rownames(MC_K) <- c("Probabilité p_n", "Erreur relative", "Itération tau")
MC_K
plot(Ks, MC_K[3,], type = "l", ylab = "Itérations nécessaires", main = "Itérations pour atteindre une erreur relative de 5%")

mc_summary <- function(x) {
    m <- mean(x)
    v <- var(x)
    e <- 1.96*sqrt(v / length(x))
    return(c("Mean" = m, "Var" = v, "Erreur abs" = e, "Erreur rel" = e / m))
}

# estimateur Monte Carlo 
K <- 22
s <- replicate(n, S())
x <- (s > K)
(r_MC <- mc_summary(x))

K <- 22
M <- function() {
    S <- 0
    M <- 0
    repeat {  #boucle infinie
        M <- M + 1
        X <- rlnorm(1, mu, sigma)
        S <- S + X
        if (S > K) return(M)  
    }
}

n <- 1e5
# estimateur avec préconditionnement 
phi <- function(m) {
    return(dpois(m, lambda) + ppois(m, lambda, lower.tail = FALSE))

}
m <- replicate(n, M())
x <- phi(m)
(r_cond <- mc_summary(x))

# estimateur Monte Carlo 
s <- replicate(n, S())
x <- (s > K)
(r_MC <- mc_summary(x))

# ratio de variance 
r_MC[2] / r_cond[2]

K <- 22
lambda <- 10
theta <- 0.7
lambda_tilde <- lambda * exp(theta)
SL <- function() {
    N <- rpois(1, lambda_tilde)
    X <- rlnorm(N, mu, sigma)
    L <- (lambda / lambda_tilde)^N * exp(lambda_tilde - lambda)
    if (sum(X) > K) return(L) 
    else return(0)
}

sl <- replicate(n, SL())
(r_IS <- mc_summary(sl))

r_MC["Var"] / r_IS["Var"] 

K <- 20
lambda0 <- 10
mu <- 0.1
sigma <- 0.3
n <- 1e5

# premier estimateur
J1 <- function(h) {
    lambda <- lambda0 + h
    s1 <- replicate(n, S())
    xph <- (s1 > K)
 
    lambda <- lambda0 - h
    s2 <- replicate(n, S())
    xmh <- (s2 > K)

    x <- (xph - xmh) / (2*h)
    return(mc_summary(x))
}
h <- c(1, 0.5, 0.1, 0.01)
res <- sapply(h, J1)
colnames(res) <- paste("h = ", h)
res

lambda <- 10
# deuxième estimateur
J2 <- function(h) {
    var_alea <- function() {
        U <- runif(1, 0, 1)
        Nmh <- qpois(U, lambda - h)
        Nph <- qpois(U, lambda + h)
        Xmh <- rlnorm(Nmh, mu, sigma)
        Xph <- rlnorm(Nph, mu, sigma)
        return(((sum(Xph) > K) - (sum(Xmh) > K)) / (2*h))
    }
    x <- replicate(n, var_alea())
    mc_summary(x)
}

h <- c(1, 0.5, 0.1, 0.01)
res <- sapply(h, J2)
colnames(res) <- paste("h = ", h)
res

lambda <- 10
# troisième estimateur
J3 <- function(h) {
    var_alea <- function() {
        U <- runif(1, 0, 1)
        Nmh <- qpois(U, lambda - h)
        Nph <- qpois(U, lambda + h)
        Nmin = min(Nmh, Nph)
        Xmin <- rlnorm(Nmin, mu, sigma)
        Xmh <- c(Xmin, rlnorm(max(Nmh-Nmin,0), mu, sigma))
        Xph <- c(Xmin, rlnorm(max(Nph-Nmin,0), mu, sigma))        
        return(((sum(Xph) > K) - (sum(Xmh) > K)) / (2*h))
    }
    x <- replicate(n, var_alea())
    mc_summary(x)
}

h <- c(1, 0.5, 0.1, 0.01, 0.001)
res <- sapply(h, J3)
colnames(res) <- paste("h = ", h)
res
