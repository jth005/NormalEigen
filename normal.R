library(CompQuadForm)
library(energy)

EYY <- function(){(2)/sqrt(pi)}

ExY <- function(dat, mu){
    n <- length(dat)
    store <- rep(0,n)
    for(i in 1:n){
       store[i] <- 2*(dat[i]) * pnorm(dat[i], mean=mu, sd=1) + 2 * dnorm(dat[i], mean=mu, sd=1) - (dat[i] - mu)
    }
    return(store)
}

norm.kernel <- function(dat, mu=NULL, sigma=NULL){
    n <- length(dat)
    if(is.null(mu) == TRUE)
        mu <- mean(dat)
    if(is.null(sigma) == TRUE){
        sigma <- sd(dat)}
    dat <- (dat - mu)/sigma
    mu <- 0
    sigma <- 1
    lilDist <- ExY(dat, mu)
    bigDist <- EYY()
    kmat <- matrix(0,n,n)
    for(i in 1:n){
        for(j in i:n){
            kmat[i,j] <- (1/n) * (lilDist[i] + lilDist[j] - bigDist - abs(dat[i] - dat[j])) 
            kmat[j,i] <- kmat[i,j]
        }
    }
    estat <- n*(2*mean(ExY(dat,mu)) - EYY() - 2/n^2 * sum((2*seq(n) - 1 - n)*sort(dat)))
    eigs <- eigen(kmat, only.values=TRUE, symmetric=TRUE)$values
    return(list(estat = estat, eigs = eigs))
}

M <- 300
save.estat <- rep(0,M)
eig.probs <- rep(0,M)
empirical.probs <- rep(0,M)

for (i in 1:M){
    dat <- rnorm(100, mean=0, sd=1)
    xx <- norm.kernel(dat, mu=NULL, sigma=NULL) 
    save.estat[i] <- xx$estat
    eig.probs[i] <- 1 - imhof(xx$estat, lambda = xx$eigs[1:length(dat)-1])$Qq
}

for (i in 1:M){
    empirical.probs[i] <- sum(save.estat <= save.estat[i])/M
}

plot(empirical.probs, eig.probs, xlim = c(0,1), ylim = c(0,1));abline(coef=c(0,1))

evnew <- c(0.11311, 0.08357, 0.03911, 0.03182, 0.01990, 
           0.01698, 0.01207, 0.01060, 0.00810,
           0.0072590,  0.0058203, 0.0052881, 
           0.0043833, 0.0040262, 0.0034204, 
           0.0031689, 0.0027437, 0.0025597, 
           0.0022498, 0.0021112, 0.0018784, 
           0.0017712, 0.0015920, 0.0015074, 
           0.0013665, 0.0012985, 0.0011857, 
           0.0011303, 0.0010386, 0.00099285)

dat <- rnorm(100, mean=1)
my.eigs<- norm.kernel(dat, mu=NULL, sigma=NULL)$eigs[1:30]
plot(my.eigs, evnew)
head(my.eigs)


m <- M
out <- replicate(m, expr={
  x <- rnorm(n)
  normal.e(x)
})

mean(out)
var(out)
