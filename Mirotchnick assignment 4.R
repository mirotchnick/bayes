library(R2jags)
setwd("/Users/Nick/Dropbox/Work/Toronto/Bayes")

#1

runif(1)
x <- c(16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46)
y <- c(2508,2518,3304,3423,3057,3190,3500,3883,3823,3646,3708,3333,3517,3241,3103,2776)

xyplot(y~x)

data <- list("x","y")
inits <- function(){list(tau=1,b1=1,b2=1)}
parameters <- c("tau","b1","b2","res","resrep","sres","yrep","psmaller","resrep","sresrep")

cat("
model{     for(i in 1:16){       y[i]~dnorm(mu[i],tau)       mu[i] <- b1 + b2*(x[i]-31)
       res[i] <- y[i] - mu[i]
       sres[i] <- res[i]*sqrt(tau)
       yrep[i]~dnorm(mu[i],tau)
       psmaller[i] <- step(y[i] - yrep[i])
       resrep[i] <- yrep[i] - mu[i]
       sresrep[i] <- resrep[i]*sqrt(tau)     }     b1~dnorm(0,.000001)     b2~dnorm(0,.000001)     tau~dgamma(.0001,.0001)     }", file="lin.txt")
     
lin.sim <- jags(data,inits,parameters,model.file="lin.txt",n.chains=3,n.iter=10000,n.burnin=700,n.thin=1)
attach(lin.sim)
attach.jags(lin.sim)
print(lin.sim)
write.csv(lin.sim$BUGSoutput$summary,"sum.csv")

min(mean(psmaller[,1]),(1-mean(psmaller[,1])))
min(mean(psmaller[,2]),(1-mean(psmaller[,2])))
min(mean(psmaller[,3]),(1-mean(psmaller[,3])))
min(mean(psmaller[,4]),(1-mean(psmaller[,4])))
min(mean(psmaller[,5]),(1-mean(psmaller[,5])))
min(mean(psmaller[,6]),(1-mean(psmaller[,6])))
min(mean(psmaller[,7]),(1-mean(psmaller[,7])))
min(mean(psmaller[,8]),(1-mean(psmaller[,8])))
min(mean(psmaller[,9]),(1-mean(psmaller[,9])))
min(mean(psmaller[,10]),(1-mean(psmaller[,10])))
min(mean(psmaller[,11]),(1-mean(psmaller[,11])))
min(mean(psmaller[,12]),(1-mean(psmaller[,12])))
min(mean(psmaller[,13]),(1-mean(psmaller[,13])))
min(mean(psmaller[,14]),(1-mean(psmaller[,14])))
min(mean(psmaller[,15]),(1-mean(psmaller[,15])))
min(mean(psmaller[,16]),(1-mean(psmaller[,16])))

like <- sum(sresrep^2)


data <- list("x","y")
inits <- function(){list(tau=1,b1=1,b2=1,b3=1)}
parameters <- c("tau","b1","b2","b3")
     
cat("
model{     for(i in 1:16){       y[i]~dnorm(mu[i],tau)       mu[i]<- b1 + b2*(x[i]-31)+ b3*pow((x[i]-31),2)     }     b1~dnorm(0,.000001)     b2~dnorm(0,.000001)     b3~dnorm(0,.01)     tau~dgamma(.0001,.0001)     }", file="quad.txt")

quad.sim <- jags(data,inits,parameters,model.file="quad.txt",n.chains=3,n.iter=10000,n.burnin=700,n.thin=1)
attach.jags(quad.sim)
attach(quad.sim)
print(quad.sim)
write.csv(BUGSoutput$summary,"sum.csv")

sx1 <- as.vector(scale(x))
sx2<- as.vector(scale(x^2))
sy <- as.vector(scale(y))


data <- list("sx1","sx2","sy")
inits <- function(){list(tau=1,b1=1,b2=1,d=1)}
parameters <- c("tau","b1","b2","d")

cat("
model{
	for(i in 1:16){
		sy[i]~dnorm(mu[i],tau)
		mu[i] <- b1*sx1[i] + d*b2*sx2[i]

	}
	d~dbern(0.5)
	b1~dnorm(0,0.25)
	b2~dnorm(0,0.25)
	tau~dgamma(.5,.01)
}", file="bf.txt")

bf.sim <- jags(data,inits,parameters,model.file="bf.txt",n.chains=3,n.iter=10000,n.burnin=700,n.thin=1)
print(bf.sim)
attach.jags(bf.sim)
mean(b1)/mean(b2)

residuals()

#2
g=function(x){(x>0)*(x<1)*((x<=0.5)*4*x+ (x>0.5)*(4-4*x))}
plot(g)

f=function(x){(x>=0)*(x<=1)}
plot(f)

#a)
u1 <- 0
u2 <- 0
z <- 0
for(i in 1:500){
	u1[i] <- runif(1)
	u2[i] <- runif(1)
	z[i] <- (u1[i]+u2[i])/2
}

u1 <- runif(500)
u2 <- runif(500)
z <- (u1+u2)/2

plot(density(z))
mean(z)
var(z)

#b)
u <- runif(1000)
w <- g(u)
plot(density(w))
1/1000*sum(w)

u <- 0
w <- 0
e <- 0
es <- 0
for(i in 1:1000){
	u[i] <- runif(1)
	w[i] <- g(u[i])
	e[i] <- u[i]*w[i]
	es[i] <- u[i]*u[i]*w[i]
}

e1 <- 1/1000*sum(e)
e2 <- 1/1000*sum(es)
v <- e2 - e1^2
var(e)

#c)
x <- runif(500)

a <- 0
for(i in 1:500){
	u <- runif(1)
	if(u <= 0.5*g(x[i]))
		a[i] <- x[i]
}

b <- na.omit(a)
length(b)/500
mean(b)
var(b)
plot(density(b))

#d)
x <- rep(0,500)
x[1] <- 0.5
ys <- 0
a <- 0
b <- 0
for(i in 2:500){
	ys[i] <- runif(1)
	tf <- g(ys[i])/g(x[i-1])
	if(runif(1) < tf){	
		x[i] <- ys[i]
		a <- a + 1
	}
	else{
		x[i] <- x[i-1]
		b <- b + 1
	}
}

a/(a + b)
mean(x)
var(x)