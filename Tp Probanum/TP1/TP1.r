x <- 17
x
y <- x + 3
y

y**2  #ou y^2

w=1:4
w

m <- matrix(w, nrow = 3, ncol = 4) 
m

m+2 

example(plot)

help(RNGkind)
RNGkind()
seed <- 12345
set.seed(seed, kind="Mersenne-Twister", normal.kind="Box-Muller")
RNGkind()
.Random.seed

help(runif)

nsim <- 10**4  #nbre de simulation
x <- runif(nsim) # vecteur de 10**4 réalisations d'une loi uniforme sur [0,1]
Xnf <- x[-nsim]       # Xn = x supprimer de sa dernière composante.
Xnd <- x[-1]        # Xnd = x supprimer de sa première composante.

par(mfrow = c(1,3))
hist(x, main="Histogramme", ylab=NA)                         #histogramme empirique 
plot(Xnf, Xnd, main="Paires adjacentes", xlab=NA, ylab=NA)   # on affiche le nuage des points (Xnf, Xnd)
acf(x, main="Fonction d'autocorrélation", ylab=NA)           # on trace la fonction d'autocorrélation

par(mfrow = c(1,2))

x <- rexp(5000)
hist(x, breaks = 50, prob = TRUE, main = "Réalisations rexp",     # histogramme normalisé (prob = T) et breaks = 50
     xlab = NA, ylab = NA, xlim = c(0, 8), ylim = c(0,1))          
curve(dexp(x), col = "red", add = TRUE)                           # on ajoute la densité de exp (add = T)

u <- runif(5000)
x <- -log(u)
hist(x, breaks = 50, prob = TRUE, main = "Réalisations -log(runif)", 
     xlab = NA, ylab = NA, xlim = c(0, 8), ylim = c(0,1))
curve(dexp(x), col = "red", add = TRUE)                                   

alpha <- 1.5
beta <- 3
c <- 2

par(mfrow=c(1,2))
plot(x = NA, y = NA, xlim = c(0,1), ylim=c(0, c), # le carré [0,1] x [0,c]  qu'on remplira
     main="Représentation du rejet", ylab=NA, xlab=NA)
curve(dbeta(x, alpha, beta), lwd=2, col="blue", add=T)      # la densité f(x) en bleu comme sur la figure
abline(h = c, lwd = 2, col="red") # cg(x) avec g(x) en rouge

y <- runif(1000)  # les (Y_n) selon la loi proposée de densité g
u <- runif(1000)   # les (U_n) pour faire le test
accept <- u*c < dbeta(y, alpha, beta)               # test du rejet 
points(y[accept], c*u[accept], pch=20, col="blue")      # points acceptés
points(y[!accept], c*u[!accept], pch=20, col="grey")  #  points rejetés

# y[accept] contient les (X_n) et doit être distribué selon la loi beta 
hist(y[accept], prob=T, xlim = c(0,1), ylim=c(0,c), main="Histogramme", xlab=NA, ylab=NA)
curve(dbeta(x, alpha, beta), lwd=2, col="blue", add=TRUE)

optim <- optimize(function(x) dbeta(x, alpha, beta),      
                interval=c(0,1), maximum=T)             
c <- optim$objective          # on recupere l'argmax dans objective 


c

help(rbinom)

rbinom_def <- function(n, size, prob) {         
    realisation <- function(size, prob) {       
        u <- runif(size)
        b <- u[u < prob]                        # u < prob est un vecteur avce des 0 et 1 
        return(length(b))
    }
    replicate(n, realisation(size, prob))  # vecteur de taille n
}

a<-rbinom_def(10,3,0.2)

a

rbinom_invF <- function(n, size, prob) {
    somme_pk <- cumsum(dbinom(0:size, size, prob))    # somme cumulée de la binomiale
    realisation <- function(somme_pk) {
        u <- runif(1)
        return(findInterval(u, somme_pk))
    }
    replicate(n, realisation(somme_pk))
}

rbinom_invF(10,3,0.2)

n <- 20  # pour n =20 sur la figure ci desssus
p <- 0.5
par(mfrow = c(1,2))
k <- 0:n
N <- 10^5

x <- rbinom_def(N, n, p)
hist_x <- table(x)/N
plot(hist_x, type="h", xlim=c(0,n), ylim=c(0,0.2), 
     xlab=NA, ylab=NA, main="rbinom_def")
points(k, dbinom(k, n, p), col = "blue", lwd = 2)

x <- rbinom_invF(N, n, p)
hist_x <- table(x)/N
plot(hist_x, type="h", xlim=c(0,n), ylim=c(0,0.2), 
     xlab=NA, ylab=NA, main="rbinom_invF")
points(k, dbinom(k, n, p), col = "blue", lwd = 2)

system.time(rbinom_def(10^5, 100, 0.5))
system.time(rbinom_invF(10^5, 100, 0.5))

nu <- c(0.1,0.3,0.3,0.3)# j'ai arrondi
n <- 1000
x <- sample(c(1:4), n, replace = T, prob = nu)
barplot(table(x)/n, ylim=c(0,0.5))  #diagramme en barre 

n <- length(nu)
P <- diag(n*nu)
S <- which(nu < 1/n)
L <- which(nu > 1/n)
I <- c()
print(P)

while (length(S) > 0) {
    # coeur de l'algorithme 
    i <- S[1]
    j <- L[1]
    delta <- 1-P[i,i]
    P[i,j] <- delta
    P[j,j] <- P[j,j] - delta
    S <- which(diag(P) < 1)
    L <- which(diag(P) > 1)
    I <- c(I, i)
    S <- setdiff(S, I)

    # on affiche
    print(S)
    print(L)
    print(P)
}

library(ggplot2)
help(ggplot)
example(ggplot)

library(cowplot)
par(mfrow = c(1,3))
df1=data.frame(x)
a<-ggplot(df1, aes(x=x)) + geom_histogram(bins=30,fill="white", color="black")+labs(title="Histogramme",x="x")
df2=data.frame(Xnf,Xnd)
b<-ggplot(df2, aes(x=Xnf, y=Xnd)) + geom_point()
conf.level <- 0.95 
ciline <- qnorm((1 - conf.level)/2)/sqrt(length(x)) 
bacf <- acf(x, plot = FALSE) 
bacfdf <- with(bacf, data.frame(lag, acf)) 

library(ggplot2) 
c<- ggplot(data=bacfdf, mapping=aes(x=lag, y=acf)) + 
     geom_bar(stat = "identity", position = "identity") 


plot_grid(a, b, c, ncol = 3, nrow = 1)

