mu <- function(x) {
    return(exp(-abs(x)/5))
}

# la fonction suivante code la simulation de Y_{n+1} à partir de Y_n = y
phi <- function(y) {
    prop <- y
    if (runif(1) < 0.5) { prop <- y-1 }
    else { prop <- y+1 }
    v <- runif(1)
    if (v < min(mu(prop) / mu(y), 1)) { return(prop) }
    else { return(y) }
}

# une trajectoire jusqu'à l'itération n
traj <- function(n) {
    result <- numeric(n)
    y <- 0
    for (i in 1:n) {
        y <- phi(y)
        result[i] <- y
    }
    return(result)
}

n <- 1000
plot(traj(n), type = "s", ylab = "Etat y_n", xlab = "Iteration n", main="Une trajectoire de la chaine de matrice Q")
Nsim <- 1000
x <- replicate(Nsim, traj(n))
barplot(table(x[n,]) / ncol(x), main=paste("Histogramme sur", Nsim, "réalisation de Y à l'itération", n))

villes <- read.table("villes.txt")
N <- nrow(villes)
factorial(N) / (2*N)

distance <- function(i, j) {
    return(sqrt((villes[i,1]-villes[j,1])^2 + (villes[i,2]-villes[j,2])^2))
}

H <- function(x) {
    sum <- distance(x[N], x[1])
    for (i in 1:(N-1)) {
        sum <- sum + distance(x[i], x[i+1])
    }
    return(sum)
}

affiche <- function(x) {
    longueur <- format(H(x), digits = 5)
    plot(villes[,1], villes[,2], col = "blue", xlab = NA, ylab = NA, main = paste("Chemin de longueur", longueur))
    points(villes[1,1], villes[1,2], col = "red")
    u <- c(villes[x, 1], villes[x[1], 1])
    v <- c(villes[x, 2], villes[x[1], 2])
    lines(u, v, col = "grey")
}

affiche(1:N)

P1 <- function(path) {
    ec <- sample(1:N, 2, replace = TRUE)
    path[ec] <- path[rev(ec)]
    return(path)
}

Q <- function(P, T, x, n) {
    y <- P(x)
    delta <- H(y) - H(x)
    if (delta < 0) { return(y) } 
    if (runif(1) < exp(-delta / T(n))) { return(y) } 
    else { return(x) }
}

T <- function(n) { return(0.01 / log(n)) }
x <- 1:N
N0 <- 1000
longueurs <- numeric(N0)
for (n in 1:N0) {
    x <- Q(P1, T, x, n)
    longueurs[n] <- H(x)
}
x
#png("chemin_P1.png", width = 800)
par(mfrow = c(1,2))
affiche(x)
plot(longueurs, type = "l", ylab=NA, xlab=NA, main="Longueur en fonction de n")

P2 <- function(path) {
    ec <- sample(1:N, 2, replace = FALSE)
    path[ec] <- path[rev(ec)]
    return(path)
}

P3 <- function(path) {
    ec <- sample(1:N, 2, replace = FALSE)
    indices <- min(ec):max(ec)
    path[indices] <- path[rev(indices)]
    return(path)
}

x <- read.table("best.txt")
affiche(x$V2)
