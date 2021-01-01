cop_gauss <- function(n, rho) {
    copule <- list(n = n, rho = rho)
    copule$rand <- function() {
        g <- sqrt(rho) * rnorm(1) + sqrt(1-rho) * rnorm(n)
        u <- pnorm(g)
        return(u)
    }     
    return(copule)
}
cop <- cop_gauss(n=2, rho = 0.)
cop$rand()    # exemple d'une réalisation du vecteur (U_1, U_2)

par(mfrow=c(1,2))
col <- rgb(0, 0, 1, alpha = 0.4)
for(rho in c(0.2, 0.8)) {
    cop <- cop_gauss(n = 2, rho = rho)
    u <- replicate(2000, cop$rand())
    plot(u[1,], u[2,], pch = 20, col = col, main = paste("rho = ", as.character(rho)))
}


randX <- function(lambda, copule) {
    u <- copule$rand()
    tau <- qexp(u, rate = lambda)            
    return(sum(tau <= 1))   #lindicatrice
}

lambda <- c(0.4)
list_rho <- c(0, 0.2, 0.6, 0.8) 
result <- matrix(nrow = 2000, ncol = 4)      # on stocke 2000 réalisations de X pour 4 valeurs de rho 
colnames(result) <- list_rho                 # on donne des noms aux colonnes pour faciliter l'accès
for (rho in list_rho) {
    cop <- cop_gauss(n = 10, rho = rho)
    x <- replicate(nrow(result), randX(lambda, cop))
    result[,as.character(rho)] <- x
}

#valeur des moyennes empiriques
apply(result, 2, mean)
print(paste("valeur théorique:", cop$n*(1-exp(-lambda))))


tables <- apply(result, 2, function(x) table(factor(x, 0:10)) / length(x)) #apply
# fonction factor(x, 0:10) x à valeurs dans 0,1,...,10. 

    
#premier graphe
par(mfrow=c(2,2))
for (i in 1:ncol(tables)) {
    barplot(tables[,i], col = rgb(0,0,1,list_rho[i]), ylim = c(0,0.4), main = paste("rho = ", list_rho[i]))
}

# second graphe
par(mfrow=c(1,1))
barplot(t(tables), col = rgb(0,0,1,list_rho), beside = TRUE, legend = paste("rho = ", list_rho), 
        main = "Evolution de la distribution de X pour différentes valeurs de rho")


Ehrenfest_dynamique <- function(state) {
    i <- sample(state$I, size = 1)
    if (i %in% state$A) { 
        state$A <- setdiff(state$A, i)
        state$B <- c(state$B, i)
    } else {
        state$B <- setdiff(state$B, i)
        state$A <- c(state$A, i)
    }
    return(state)
}

# on simule 1 realisation de la chaine de Markov
Ehrenfest_realisation <- function(state, n) {
    for (k in 1:n) {
        state <- Ehrenfest_dynamique(state)
    }
    return(length(state$A))  #on renvoie X_n
}

init_state <- list(I = 1:40, A = 1:10, B = 11:40)
x <- replicate(5000, Ehrenfest_realisation(init_state, 100))
barplot(table(x) / length(x))

x_0 <- 10
Ehrenfest_transition <- function(x_n) {
    if (x_n == 0) return(1)
    if (x_n == d) return(d-1)                   
    delta <- sample(c(-1, +1), 1, prob = c(x_n/d, (d-x_n)/d))
    return(x_n + delta) #x_n+1
}

Ehrenfest_realisation <- function(x, n) {
    for (k in 1:n) {
        x <- Ehrenfest_transition(x)
    }
    return(x)
}

par(mfrow=c(1,2))
d <- 40
x <- replicate(1000, Ehrenfest_realisation(x_0, 100))    #pair 
barplot(table(x) / length(x), main="Histogramme de X_100")
x <- replicate(1000, Ehrenfest_realisation(x_0, 101))    #impair 
barplot(table(x) / length(x), main="Histogramme de X_101")



Npaths <- 200
lambda <- 20
T <- 1
plot(NULL, xlim=c(0,T), ylim=c(0,2*lambda), xlab="Time t", ylab="Processus de comptage N_t",
     main = paste(Npaths, "trajectoires d'un processus de Poisson de paramètre", lambda))

# boucle for et d'une boucle while
for (k in 1:Npaths) {
    T_n <- 0
    path <- 0
    while (T_n < T) {
        T_n <- T_n + rexp(1, lambda)
        if (T_n < T) path <- c(path, T_n)
    }
    N_T <- length(path)-1
    N <- 0:N_T

    # on ajoute le point (T, N_T)    
    path <- c(path, T)
    N <- c(N, N_T)
    
    lines(path, N, type="s", col=rgb(0,0,1,alpha=0.5))
}


Npaths <- 200
lambda <- 20
T <- 1
plot(NULL, xlim=c(0,T), ylim=c(0,2*lambda), xlab="Time t", ylab="Processus de comptage N_t",
     main = paste(Npaths, "trajectoires d'un processus de Poisson de paramètre", lambda))

for (k in 1:Npaths) {
    N_T <- rpois(1, lambda*T)
    T_n <- runif(N_T, min = 0, max = T)
    T_n <- c(0, sort(T_n), T)
    lines(T_n, c(0:N_T, N_T), type="s", col=rgb(0,0,1,alpha=0.5))
}

#code plus court 




