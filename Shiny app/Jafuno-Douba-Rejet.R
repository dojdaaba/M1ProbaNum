#' ---	
#' title: "Méthode de Rejet"	
#' author: "Douba Jafuno"	
#' date: "5 avril 2019"	
#' output: html_document	
#' runtime: shiny	
#' ---	
#' 	
#' 	
knitr::opts_chunk$set(echo = TRUE)	
#' 	
#' 	
#' ## Introduction 	
#' 	
#' ## But	
#' 	
#' Que faisons-nous si nous voulons générer des échantillons d'une variable aléatoire X de densité $f(x)$, que nous ne sachions pas simuler directement ? 	
#' 	
#' Avec la **Methode de rejet**  bien que nous ne puissons pas facilement échantillonner $f(x)$, il éxiste une autre densité $g(x)$, à partir de laquelle il est plus facile pour nous d'échantillonner.	
#' 	
#' Nous pouvons ensuite échantillonner directement dans $g$, puis rejeter les échantillons de manière stratégique pour donner l'impression que les échantillons non rejetés résultants ressemblent à la densité $f$	
#' 	
#' La densité $g$ sera appellé densité candidate et $f$ sera la densité cible.	
#' 	
#' Soit donc $(Y,U)$ un couple de variables aléatoires indépendantes  tirées selon une loi uniforme avec $Y\sim g$ , i.e. $(Y,U)$  est un point tiré uniformément dans le carré unité par exemple. On peut alors montrer que la distribution de $X$ est la loi conditionnelle de $Y$ sachant l'évènement $M= \{U\ge f_{X}(Y)\}$. Autrement dit, ${f_{X}(x)=f_{Y}(x|M)}$.	
#' Pour simuler une suite de variables aléatoires réelles $(X_{n})_{n\ge 1}$ de distribution identique à celle de $X$ il suffit donc, dans une suite de tirages de couples $(Y_{i},U_{i})$ uniformes indépendants, de sélectionner les $Y_{i}$ correspondant aux tirages $(Y_{i},U_{i})$ vérifiant $M$, et de rejeter les autres. Tout celà donne lieu à un algorithme.	
#' 	
#' ## Algorithme	
#' Pour utiliser l'algorithme de la méthode des rejets, nous devons d'abord nous assurer que le support de $f$ est un sous-ensemble du support de $g$	
#' . Si  $X_{f}$ est le support de $f$ et $X_{g}$ est le support de $g$, alors nous devons avoir $X_{f}\subset X_{g}$	
#' 	
#' Cela a du sens: s'il y a une région du support de $f$ que $g$ ne peut jamais toucher, alors cette région ne sera jamais échantillonnée. De plus, il faut supposer que $c$ = $\sup_{x\subseteq X_{f}}\frac{f(x)}{g(x)} < \infty$ 	
#' 	
#' L'algorithme de la Méthode de rejet pour prélever un échantillon à partir de la densité cible $f$ est alors: 	
#' 	
#' - Simuler $U\sim \mathcal{U}([0,1])$	
#'   	
#' - Simuler $Y\sim g$	
#'   	
#' - Si $U < \frac{f(Y)}{cg(Y)}$ 	
#'   On accepte Y sinon on rejette Y et on recommence depuis le début.	
#' 	
#' De plus le nombre N d'itération de l'algorithme suit une loi géométrique de paramètre $\frac{1}{c}$ 	
#' Ce paramètre correspond donc à la probabilité d'accepter Y en effet:	
#' 	
#' $P(Y:accept?)=P(U < \frac{f(Y)}{cg(Y)})=\int_{}P(U < \frac{f(y)}{cg(y)})|Y=y)g(y)dy= \int_{}\frac{f(y)}{cg(y)}g(y)dy=\frac{1}{c}$	
#' 	
#' ## Exemple: Methode de rejet pour la loi gamma	
#' 	
#' Soit  X $\sim$ $\Gamma(\alpha,\beta)$, Loi Gamma de paramètre $\alpha$ > 0 et $\beta$ > 0 de densité $f_{\alpha, \beta}(x) = \frac{1}{\Gamma(\alpha)} \beta^{\alpha}x^{\alpha-1} e^{-\beta x} \mathbf{1}_{]0,+\infty[}(x)$	
#' 	
#' Comment Simuler $X$?	
#' 	
#' - Si $\alpha$ = 1, cela coincide avec une variable aléatoire de loi exponnentielle de parametre $\alpha$.	
#' 	
#' - Si $\alpha$ est un entier non nul, une variable aléatoire de loi $\Gamma(\alpha,\beta)$ s'écrit comme la somme de $\alpha$ variables aléatoires indépendantes de loi exponnentielle de paramètre $\beta$ : la simulation en découle immédiatement (Méthode d'inversion)	
#' 	
#' - Si $\alpha$ n'est pas entier, il est utile d'avoir recours à la **Méthode de Rejet**. Donnons en une illustration, sans avoir le soucis d'optimalité. Pour simplifier, supposons $\beta$ = 1 et $\alpha \in (n, n + 1)$ 	
#' 	
#' Pour le dernier cas $\alpha$ n'est pas un entier tel que $\alpha \in (n, n + 1)$  et $\beta =1$ 	
#' 	
#' - Si $\alpha > 1$. Prenons pour $Y \sim g$  une loi **$\Gamma(n,\frac{1}{2})$**	
#' .On vérifie alors que la constante de rejet vaut:  $c = \sup_{x>0}\frac{f(x)}{g(x)} =\frac{\Gamma(n)2^{\alpha}(\alpha-n)^{\alpha-n}e^{-(\alpha-n)}}{\Gamma(\alpha)}$.Cette constante augmente rapidement lorsque $\alpha$ tend vers +$\infty$.	
#' 	
#' 	
#' 	
#' - Si $\alpha \in ]0,1[$ et donc $\beta=1$ avec $f_{\alpha, 1}(x) = \frac{1}{\Gamma(\alpha)}x^{\alpha-1} e^{-x} \mathbf{1}_{]0,\infty[}(x)$. Prenons pour $Y \sim g$ de loi que l'on surnomme **Auxgamma** de densité $g(x)$ = $\frac{\alpha e}{\alpha +e}(x^{\alpha-1}\mathbf{1}_{]0,1[}(x) + e^{-x}\mathbf{1}_{[1,+\infty[}(x))$. Pour simuler selon la loi Auxgamma de densité g on utilise l'inverse de la fonction de repartition. La fonction quantile s'écrit $\forall U \in ]0,1[$,  $G^{-1}(u)=(\frac{(\alpha +e)u}{e})^{\frac{1}{\alpha}}\mathbf{1}_{u<\frac{e}{\alpha +e}} - log((1-u)\frac{(\alpha +e)u}{e})\mathbf{1}_{u\ge\frac{e}{\alpha +e}}$ 	
#' .On vérifie alors que la constante de rejet vaut:  $c = \frac{(\alpha +e)}{\alpha e\Gamma(\alpha)}$	
#' 	
#' ## Application 	
#' 	
shinyApp(	
	
	
ui <- fluidPage(	
  titlePanel("Methode de Rejet"),	
  	
  sidebarLayout(	
   sidebarPanel(	
     sliderInput("n","Taille de l'echantillion:",min = 0,max = 10e3,value = 1000, step = 50),	
    textInput("f", "Entrez la densite f(x) a simuler 	
                (vous pouvez utiliser la documentation de R):", "dgamma(x,0.5,1)"),	
      tags$h6("Par exemple entrer dnorm(x,0,1) pour une loi normale centree reduite"),	
      verbatimTextOutput("value"),	
      numericInput("c", "Entrer la constante C (affiche dans la deuxieme table) ", 1.4),	
      tags$b("Choissez une loi de densite g(x) parmi:"),	
     	
      	
     div(  tabsetPanel(id = "tabset",	
                  	
                  tabPanel("AuxGamma",	
                           tags$h6("Utilisez cette loi pour simuler une Gamma de parametre Alpha et 1 avec Alpha dans ]0,1[" ),	
                           numericInput("Alpha", "Alpha", 0.5)	
                           ),	
                  tabPanel("Gamma",	
                           tags$h6("Vous pouvez l'utiliser pour simuler une Gamma(alpha,beta) avec Alpha > 1 ici entrer abs(Alpha) et 0.5"),	
                           numericInput("alpha", "alpha", 1),	
                           numericInput("beta", "beta", 0.5)	
                  ), 	
                  tabPanel("Weibull",	
                           	
                           numericInput("a", "a", 0.5)	
                          	
                          	
                  ),	
                  tabPanel("Uniforme",	
                           	
                           sliderInput("uni_range", "Loi uniforme sur [a,b] ", min = -10, max = 10, value = c(0, 1))	
                  ),	
                          	
                  	
                 	
                  tabPanel("Gaussienne",	
                           	
                           numericInput("mu", "Mu", 0),	
                           numericInput("sigma", "Sigma", 1)	
                  )	
      ),class = "span1"),	
      	
          actionButton("go", "Demarrer"), width ='4'),	
  	
    mainPanel(	
      tabsetPanel(id = "tabset",	
            tabPanel("Representation et Histogramme", plotOutput("plot")), 	
  sliderInput("plotrange", "axe x", min = -10, max = 10, value = c(0, 4), step = 0.1)	
     	
      )	
    ))),	
	
	
server <- function(input, output) {	
  observeEvent(input$go,{	
    	
    rwei<-function(n,a){	
      wei<-function(a){	
        u<-runif(1)	
        G=(-log(1-u))^(1/a)	
        return(G)	
      }	
      return(replicate(n,wei(a)))	
    }	
    	
    dwei<-function(x,a){	
      return (x^(a-1)*exp(-x^a))	
    }	
    rauxgamma<-function(n,alpha){	
      auxgamma<-function(alpha){	
        u<-runif(1)	
        g<-(alpha +exp(1))/exp(1)	
        G<-(g*u)^(1/alpha)*(u<(1/g))-log((1-u)*(g/alpha))*(u>=(1/g))	
        return(G)	
        	
      }	
      return (replicate(n, auxgamma(alpha)))	
      	
    }	
    	
    dauxgamma<-function(x,alpha){	
      g<-(alpha*exp(1))/(alpha +exp(1))	
      G<- g*(x^(alpha-1)*dunif(x)+exp(-x)*(x>=1))	
      return (G)	
    }	
    	
    	
    	
    	
    	
    accept_reject1 <- function(n) {	
     	
        if (input$tabset == "Gaussienne") {	
          x <- rnorm(n,input$mu, input$sigma)       	
          U <- runif(n, 0, 1)                   	
          y <- input$c*U*dnorm(x,input$mu, input$sigma)	
        }	
        if (input$tabset == "Uniforme") {	
          x <- runif(n, input$uni_range[1],input$uni_range[2])       	
          U <- runif(n,  0,  1)*dunif(x, input$uni_range[1], input$uni_range[2])    	
          y <- input$c*U	
        }	
        if (input$tabset == "AuxGamma") {	
          x <- rauxgamma(n, input$Alpha)          	
          U <- runif(n, 0, 1)*dauxgamma(x, input$Alpha)	
          y <- input$c*U	
        }	
        if (input$tabset == "Gamma") {	
          x <- rgamma(n, input$alpha,  input$beta)           	
          U <- runif(n, 0, 1)*dgamma(x, input$alpha, input$beta)	
          y <- input$c*U	
        }	
       if (input$tabset == "Weibull") {	
        x <- rwei(n, input$a)           	
        U <- runif(n, 0, 1)*dwei(x, input$a)	
        y <- input$c*U	
        }	
        accept<- y <= f(x) 	
        	
        a <- x[accept] # x[c*u*g(x)<f(x)]	
        b<- y[accept] # c*u*g(x)[c*u*g(x)<f(x)]	
        c<-x[!accept] 	
        d<-y[!accept]	
      	
      	
      return(list(accept = cbind(a, b),	
                  reject = cbind(c, d)))	
    }	
     	
    	
    	
    eval(parse(text = paste('f <<- function(x)', input$f, sep='')))	
    if (input$tabset == "Gaussienne") {	
      op<- optimize(function(x) f(x)/dnorm(x,input$mu, input$sigm) ,      	
                    interval=input$plotrange, maximum=T)$objective  	
    }	
    if (input$tabset == "Uniforme") {	
      	
      op<- optimize(function(x) f(x)/dunif(x, input$uni_range[1], input$uni_range[2]) ,      	
                    interval=input$plotrange, maximum=T)$objective	
      	
    }	
    if (input$tabset == "Gamma") {	
      	
      op<- optimize(function(x) f(x)/dgamma(x, input$alpha, input$beta) ,      	
                    interval=input$plotrange, maximum=T)$objective	
      	
    }	
    	
    	
    	
    if (input$tabset == "AuxGamma") {	
      	
      op<- optimize(function(x) f(x)/dauxgamma(x, input$Alpha) ,      	
                    interval=input$plotrange, maximum=T)$objective	
      	
    }	
    	
    if (input$tabset == "Weibull") {	
      	
      op<- optimize(function(x) f(x)/dwei(x, input$a) ,      	
                    interval=input$plotrange, maximum=T)$objective	
      	
    }	
    	
    l <- accept_reject1(input$n)	
    accept <- l$accept	
    reject <- l$reject	
   	
    output$plot <- renderPlot({	
      par(mfrow=c(1,2))	
      plot(NA, xlim = input$plotrange,	
           ylim = c(0,input$c), main="Representation du Rejet ",xlab = paste("Constante c = ", op), ylab = "",axes = TRUE)	
      curve(f(x), lwd=2, col="blue", add=TRUE)	
      	
      abline(h = input$c, lwd = 2, col="red") 	
      points(accept[,1],accept[,2],pch= 20, col = "blue") 	
      points(reject[-1,1], reject[-1,2], pch=20, col = "grey")	
      legend(input$plotrange[2]-input$plotrange[2]*30/100, input$c, c('accepte', 'rejete'),	
             col = c('blue', 'grey'), pch = c(20, 20))	
      legend(input$plotrange[2]-input$plotrange[2]*40/100, input$c-input$c*20/100, c('Constante c'),lty=1,	
             col = 'red')	
     	
      hist(accept[,1], prob=T,xlim = input$plotrange,ylim = c(0,input$c), main="Histogramme", xlab=NA, ylab=NA)	
      curve(f(x), lwd=2, col="blue", add=TRUE)	
      	
    })	
    	
  })	
},	
options = list(	
    width = 1000, height = 900)	
)	
	
#' 	
#' 	
