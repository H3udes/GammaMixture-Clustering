################################################################################
#
#       Title:  Gamma mixtures models Library
#       Autor:  Heudes Emmanuel Ochoa Parra
#       Date: 12/12/2016
#       Desc: A version of EM and SEM algorithm to fit parameters of a
#       gamma mixture. The main idea is to estimate the paremeters of
#       the associated gamma mixture model.#
#       
#       Instructions to run:
#         - Put all the data sets file in the working directory with the
#           "functionszc.R" inside.
#         - Install the needed packages
#         - Execute the code
#
################################################################################

require(mclust)
require(plotly)
#EM & SEM GAMM Function:-------------------
GammaClust <- function(data,k,Eps = 1e-4, maxiter = 1e3, conv = FALSE, beta1 = FALSE){
    
    # ----| Help: GammaClust |----
    
    #Function taht uses the "Stocastic Expectation Maximization Algorithm" (SEM)
    #to find the parameters of a "k" gamma mixture model
    
    # Arguments:
    
    # - data: data set to clusterize (Vector). 
    # - k: Number of clusters (densities) in the mix.
    # - Eps: Covengence treshold
    # - maxiter: max interations
    # - conv: If TRUE, Offers outputs toobserve the convergence
    # - beta1: If TRUE, set the beta parameter o f gamma densities in 1.
    
    # ----| Processing |----
    
    if (conv == TRUE){
        
        a <- list()
        b <- list()
        deltal <- list()
        ww <- list()
    }
    
    #   Initialize the weights, and alpha, beta parameters
    
    y <- data
    
    # Optimize the initialization of aramters using the results of a gaussian mixture model
    # that estimates the mean and variance of the data set.
    
    clust0 <- Mclust(data = y,G = k, modelNames = 'V', prior = priorControl(scale = c(sd(y))))
    
    Ey <- c()
    for (i in 1:k){
        nk <- sum(clust0$classification == i)
        Ey[i] <- mean(y[clust0$classification == i])
    }
    pos <- order(Ey,decreasing = TRUE)
    Ey <- Ey[pos]
    
    w <- clust0$parameters$pro[pos]
    
    alpha <- Ey
    beta <- rep(1,k)
    
    condS<-NA
    m <- 1 
    if ((sum(is.nan(alpha)) == 0) & (sum(is.nan(beta)) == 0) & (sum(is.infinite(alpha)) == 0) & (sum(is.infinite(beta)) == 0) & (sum(alpha < 0) == 0) & (sum(beta < 0) == 0)){
        
        # Initializing the likehood function:
        
        Lm <- eval.logverosim.gamma(y,w,alpha,k,beta)
        Lmm <- 0.8
        delta <- abs(Lmm - Lm)
        while((Eps < delta)&(m < maxiter)){
            
            #E Step
            
            E.step <- gamma.step.E(y,w,alpha,beta,k)
            
            #S Step
            
            S.step <- gamma.step.S(k,E.step$p.gamma)
            condS <- S.step$cond
            if (!condS){
                break
            }
            
            # M Step
            
            M.step <- gamma.step.M(y,k,p.gamma = S.step$p.gamma,n.j = S.step$n.j,alpha,beta)
            
            w <- as.numeric(M.step$w.n)
            alpha <- as.numeric(M.step$alpha.n)
            beta <- as.numeric(M.step$beta.n)
            
            if ((sum(is.nan(w)) > 0 ) | (sum(is.nan(alpha)) > 0) | (sum(is.nan(beta)) > 0) | (sum(is.infinite(alpha)) > 0) | (sum(is.infinite(beta)) > 0) | (sum(alpha < 0) > 0) | (sum(beta < 0) > 0)){
                break 
            }
            
            Lmm <- eval.logverosim.gamma(y,w,alpha,k,beta)
            delta <- abs(Lmm - Lm)
            Lm <- Lmm
            
            if (conv == TRUE){
                
                a[[m]] <- alpha
                b[[m]] <- beta
                deltal[[m]] <- delta
                ww[[m]] <- w
            }
            m <- m+1
        }
    }
    
    
    mu <- alpha/beta
    pos <- order(mu,decreasing = TRUE)
    mu <- mu[pos]
    
    alpha <- alpha[pos]
    beta <- beta[pos]
    
    
    
    sigmasq <- alpha/beta^2
    modes <- (alpha - 1)/beta
    if (!(sum(is.nan(alpha)) > 0 | sum(is.nan(beta)) > 0) & !(sum(is.infinite(alpha)) > 0 | sum(is.infinite(beta)) > 0)){
        if(sum(alpha < 1)){
            modes[which(alpha < 1)] <- NaN
        }
        
    }else{
        alpha[which(is.infinite(alpha))] <- NaN
        beta[which(is.infinite(beta))] <- NaN
        E.step <- gamma.step.E(y,w,alpha,beta,k)
        condS <- NA
    }
    
    
    p.gamma <- E.step$p.gamma
    p.gamma <- p.gamma[pos]
    
    #Defining the clusters according to the mixture
    
    clustmix <- c(rep(2,nrow(p.gamma)))
    pos.c1 <- p.gamma$V1 > p.gamma$V2
    pos.c2 <- p.gamma$V1 < p.gamma$V2
    
    clustmix[pos.c2] <- 1
    
    n1 <- sum(clustmix == 1)
    n2 <- sum(clustmix == 2)
    N <- n1 + n2 
    
    # Computing the uncertainity
    
    if (sum(is.na(pos.c1))>0){
        uncertainty <- rep(x = NA,nrow(p.gamma))
    }else{
        
        uncertainty <- c(1:nrow(p.gamma))
        uncertainty[pos.c2] <- p.gamma$V1[pos.c2]
        uncertainty[pos.c1] <- p.gamma$V2[pos.c1] 
    }
    
    
    # Output
    parameters <- list(w,alpha,beta,mu,sigmasq,modes,N,m-1,condS)
    
    names(parameters) <- c('w','alpha','beta','mu','sigmasq','modes','N','Iter','cond')
    
    clust <- list(clustmix,n1,n2)
    names(clust) <- c('clustmix','n1','n2')
    
    output <- list(parameters,clust,p.gamma,uncertainty)
    names(output) <- c('parameters','clust','p.gamma','uncertainty')
    
    if (conv == TRUE){
        
        a <- plyr::ldply(a,rbind)
        b <- plyr::ldply(b,rbind)
        delta <- plyr::ldply(deltal,rbind)
        ww <- plyr::ldply(ww,rbind)
        
        conv <- list(a,b,ww,delta)
        names(conv) <- c('Alpha.conv','Beta.conv','W.conv','Delta.LL')
        
        output <- list(parameters,clust,p.gamma,uncertainty, conv)
        names(output) <- c('parameters','clust','p.gamma','uncertainty','list.conv')
        return(output)
    }
    
    return(output)
    
}

# Newton's algorithm to solve the equations inside of the gamma mixture model

newton.gamma <- function(c,d,loga = FALSE,tol=1E-12,N=100) {
    
    
    # ----| Help: newton.gamma |----
    
    # Function that use the newtons algorithm to find the roots of the digamma
    # funtion.
    
    # Arguments:
    
    # - c: constant of the digamma function.
    # - d: root dimension
    # - tol: Convergence treshold
    # - N: Max iteration number
    # - loga: If TRUE it use the logarithm version.
    
    # ----| Procesamiento |----
    
    h <- rep(1e-7,d)
    i <- 1
    x0 <- rep(1.40,d)
    x1 <- x0
    p <- list()
    while (i<=N) {
        df.dx <- (digamma.c(x0+h,c,loga = loga)-digamma.c(x0,c,loga = loga))/h
        x1 <- (x0 - (digamma.c(x0,c,loga = loga)/df.dx))
        p[[i]] <- x1
        i <- i + 1
        if(sum(is.nan(x1)) > 0) break
        if (norm_vec(x1-x0) < tol & sum(x1 > 0) == d) break
        x0 <- x1 
    }
    return(x1)
}  

norm_vec <- function(x){sqrt(sum(x^2))}

# Evaluates the log-likelihood functions:------------------------

eval.logverosim.gamma <- function(y,w,alpha,k,beta){
    
    #Function that evaluates the log-likelihood function
    
    #Arguments: 
    
    # y: data
    # k: number of clusters
    # w: weights
    # alpha: shape parameter
    # beta: rate parameter
    
    #Ouput:
    
    #L: Likelihood value
    
    n <- length(y)
    w.phi <- list()
    
    for (j in 1:k ){
        
        # Evaluating the densities and weights according to alpha and beta values
        
        w.phi[[j]] <- t(w[j] * dgamma(x = y,shape = alpha[j], rate = beta[j]))
        
    }
    
    w.phi <- plyr::ldply(w.phi,rbind)
    
    w.phi <- colSums(w.phi) #Internal sumatories of the densisties and evaluated weights for each observation
    w.phi <- sapply(w.phi,log) #Computing the logarithm of each last sumatories
    L <- sum(w.phi)/n #Computing the log-likelihood
    
    return(L)
    
}


#E Step, for GAMM (Gammas Mixture):--------------------------------------

#Functiong that computes the probability of the i observation belongs to the j densitie
#Arguments: 

# y: data
# k: number of clusters (densities)
# w: weights
# alpha: shape parameter
# beta: rate parameter

#Output:

#p.gamma: Probability of the i observation belongs to the j densitie
#n.j: Numbers of elementes in the densitie j.
gamma.step.E <- function(y,w,alpha,beta,k){
    

    w.phi <- list()
    for (j in 1:k ){
        
        # Evaluating the densities and weights according to alpha and beta values
        
        w.phi[[j]] <- t(w[j] * dgamma(x = y,shape = alpha[j], rate = beta[j]))
        # w.phi[[j]] <- w.phip[[j]]
    }
    
    w.phiTotal <- colSums(plyr::ldply(w.phi,rbind)) #Total probability
    
    p.gamma <- list()
    for (j in 1:k ){
        
        # Evaluating the densities and weights according to alpha and beta values
        
        p.gamma[[j]] <- w.phi[[j]]/w.phiTotal
        
    }
    p.gammap <<- as.matrix(plyr::ldply(p.gamma,rbind)) #Probability matrix of gamma densities
    p.gamma<- p.gammap 
    
    n.jp <<- colSums(t(p.gamma))
    n.j <- n.jp
    
    p.gamma <- as.data.frame(t(p.gamma))
    n.j <- as.data.frame(t(as.matrix(n.j)))
    
    
    output <- list(p.gamma,n.j)
    names(output)<- c('p.gamma','n.j')
    
    return(output)
    
}

#S Step, for GAMM (Gammas Mixture):-----------------

#Functiong that computes the probability of the i observation belongs to the j densitie
# Using a stocastic filter trough a multinomial distribution, generating markov chains in the process
#Output: 

# p.gamma: data
# k: Number of clusters (densities)

#Output:

#p.gamma: Probability of the i observation belongs to the j densitie (filtered)
#n.j: Numbers of elementes in the densitie j.

gamma.step.S <- function(k,p.gamma,d = 1){
    cond <- TRUE
    zij <- list()
    N <- nrow(p.gamma)
    
    for (j in 1:N){
        zij[[j]] <- t(rmultinom(1,1,prob = p.gamma[j,]))
    }
    
    zij <- plyr::ldply(zij,rbind)
    
    Ccheck <- c()
    Cn <- (d+1)/N
    for (j in 1:k){
        Ccheck[j] <- sum(zij[,j])/N
    }
    
    if (sum(Ccheck >= Cn) != k){
        cond <- FALSE
    }
    
    n.jp <- colSums((zij))
    n.j <- as.data.frame(t(as.matrix(n.jp)))
    
    p.gamma <- zij
    output <- list(p.gamma,n.j,cond)
    names(output)<- c('p.gamma','n.j','cond')
    return(output)
}

#M Step, for GAMM (Gammas Mixture):------------------

#Function that computes new weight, alpha and beta parameters on each iteration.
#Arguments: 

# y: data
# k: Number of clusters (densities)
#p.gamma: Probability of the i observation belongs to the j densitie (filtered)
#n.j: Numbers of elementes in the densitie j.
# alpha: shape parameter
# beta: rate parameter


#Output:

#w.n: new weights
#alpha.n: new alpha
#beta.n: new beta

gamma.step.M <- function(y,k,p.gamma,n.j,alpha,beta){
    
    
    n <- length(y)
    
    # Computing new weights
    
    w.n <- n.j/n
    
    #Computing new alpha:
    
    suma.plogy.nj <- (colSums(p.gamma * log(y)))/n.j 
    
    
    Const <- as.numeric(log(beta) + suma.plogy.nj)
    alpha.n <- newton.gamma(c = Const,d = length(Const))
    
    #Computing new beta:
    
    alpha.nj <- alpha.n * n.j
    suma.py <- colSums(p.gamma * y) 
    
    beta.n <- alpha.nj/suma.py
    
    output <- list(w.n,alpha.n,beta.n)
    names(output)<- c('w.n','alpha.n','beta.n')
    
    return(output)
}

#Digamma functions used to find the alpha parameter:-----------------
digamma.c<- function(x,c,loga = FALSE){
    
    # ----| Help: digamma.c |----
    
    #Function that evaluates the digamma function. If it is nescesary, the function
    # has the option of evaluate a variant of the digamma function with a constant:
    #           y = digamma(x) - log(x) - c
    #Where "c" is a constant.
    
    
    # Arguments:
    
    # - x: Value to evaluate in the digamma function.
    # - c: It is the constant to use.
    # - loga: If TRUE it use the logarithm version.
    
    # ----| Processing |----
    
    k <- length(x)
    if (sum(x > 0) == k){
        if (loga == TRUE){
            
            y <- digamma(x) - log(x) - c
            
        }else{
            y <- digamma(x) - c
        }
    }else{
        y <- rep(NaN,k)
    }
    return(y) 
}

#Gamma mixture generator for two densities:-------
generator.gamma.mixturek2 <- function(a0,b0,w0,n){
    
    # ----| Help: generator.gamma.mixturek2 |----
    
    #Function that generates gamma mixtures of two densities
    
    # Arguments:
    
    # - a0: shape parameter vector
    # - b0: rate parameter vector
    # - w0: weight parameter vector
    # - n: Number of elements 
    
    # ----| Processing |----
    y <- c() 
    D1 <- 0
    D2 <- 0
    for (i in 1:n){
        piv <- runif(1)
        if(piv<=w0[1]){
            y[i] <- rgamma(1,shape = a0[1],rate = b0[1])
            D1 <- D1 + 1
        }else{
            y[i] <- rgamma(1,shape = a0[2],rate = b0[2])  
            D2 <- D2 + 1
        }
    }
    output <- list(y,c(D1,D2))
    names(output) <- c('y','D1.D2')
    return(output)
}


