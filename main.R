################################################################################
#
#       Title:  Gamma Clustering - Using Gamma Mixture Models 
#               The EM and SEM algorithm
#       Autor:  Heudes Emmanuel Ochoa Parra
#       Date: 12/12/2016
#       Desc: A version of EM and SEM algorithm to fit parameters of a
#       gamma mixture. The main idea is to estimate the paremeters of
#       the associated gamma mixture model.
#       
#       Instructions to run:
#         - Install the needed packages
#         - Execute the code
#
################################################################################

# Install required packages:
    install.packages("mclust")
    install.packages("plotly")

# Sourcing lib file:
    
    source("libgammaclust.R")
    
# Let's generate some data to fit:
    
    alpha <- c(50,60)                        #Setting alpha parameters of sample
    beta <- c(1,1)                         #Setting beta parameters of sample
    w1 <- 0.30                          #Setting weight parameters of sample
    w2 <- 1 - w1
    n <- 1000                           #Let's generate about 1000 of values
    data <- generator.gamma.mixturek2(a0 = alpha,b0 = beta,w0 = c(w1,w2),n = n)
    aux <- seq(alpha[1],alpha[2],length.out = 1000)
    data <- data.frame(x=aux,Data.Values=data$y)
    
# Let's take a look of the sample data.
    
    p <- plot_ly(data = data, x = ~x, y = ~Data.Values,type = 'scatter' ,mode = 'markers')
    p
    
# Now, let's run the gamma clustering algorithm:
    
    clusters <- GammaClust(data = data$y,k = 2)
    
    names <- c('w','alpha','beta','mu','sigmasq','modes')
    
    param.summary<-clusters$parameters[[1]]
    for (i in 2:6){
        param.summary <- cbind(param.summary,clusters$parameters[[i]])
    }
    colnames(param.summary) <- names
    View(param.summary)
    
# Looking the clustering:
    clust <- data.frame(x= aux, Data.Values=data$y,cluster=as.character(clusters$clust$clustmix))
    clust$cluster <- as.character(clust$cluster)
    clust[clust$cluster == 1,3] <- c('Lowest Density')
    clust[clust$cluster == 2,3] <- c('Highest Density')
    pal <- c('blue','red')
    c.p <- plot_ly(data = clust,x = ~aux,y = ~Data.Values,type = 'scatter' ,mode = 'markers',
                   color = ~cluster, colors = pal)
    c.p
    
    p <- plot_ly(data = iris, x = ~Sepal.Length, y = ~Petal.Length, color = ~Species)
    
    
