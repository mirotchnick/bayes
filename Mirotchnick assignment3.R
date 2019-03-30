#1

library(R2jags)
setwd("/Users/Nick/Dropbox/Work/Toronto/Bayes")

null=runif(1)
ratdat<- read.table("BigRatDat.txt")
attach(ratdat)
N<-length(V1)
wgt<-ratdat[,3:13]
dose<-ratdat[,1]
xbar<-mean(V1)

data<-list("N","wgt","dose")
inits<-function(){ list(b0=runif(50),b1=runif(50),sau=1,sb0=runif(1,0.5,1),sb1=runif(1,0.5,1))}
parameters<-c("tau","mb00","mb10","mb0diff","mb1diff","sau","sb0","sb1","tb0","tb1","b0","b1")

cat("
model{
	for(i in 1:N){   		for(j in 1:11){     		wgt[i,j]~dnorm(mu[i,j],tau)     		mu[i,j]<-b0[i]+b1[i]*j		}	b0[i]~dnorm(mb0[i],tb0)	b1[i]~dnorm(mb1[i],tb1)	mb0[i]<-mb00+mb0diff*dose[i]	mb1[i]<-mb10+mb1diff*dose[i]
	}
mb00~dnorm(100,.00001)mb10~dnorm(0,.0001)mb0diff~dnorm(0,.0001)mb1diff~dnorm(0,.0001)
tau<-pow(sau,-2)
tb0<-pow(sb0,-2)
tb1<-pow(sb1,-2)sau~dunif(0,250)sb0~dunif(0,250)sb1~dunif(0,250)
}", file="weight.txt")

weights.sim<-jags(data,inits,parameters,model.file="weight.txt", n.chains=3, n.iter=10000, n.burnin=700, n.thin=1)
attach(weights.sim)
print(weights.sim)

#b)
traceplot(weights.sim)
autocorr.plot(weights.sim,ask=T)

plot(BUGSoutput$sims.list$mb0diff[27400:27900],type="l",main="mb0diff")
plot(BUGSoutput$sims.list$mb1diff[27400:27900],type="l",,main="mb1diff")
plot(BUGSoutput$sims.list$tb0[27400:27900],type="l",main="tb0")
plot(BUGSoutput$sims.list$tb1[27400:27900],type="l",main="tb1")






write.csv(BUGSoutput$summary,"summary.csv")

#d)
plot(density(BUGSoutput$sims.list$mb0diff))
plot(density(BUGSoutput$sims.list$mb1diff))



#2
smokedie<-read.csv("smoke.csv",header=T)
attach(smokedie)
a<-length(age)

data<-list("a","age","smoke","deaths","years")
inits<-function(){list(ba=rnorm(5),bs=rnorm(4,0,1),std=runif(1,0.5,2),stda=runif(1,0.5,2),stds=runif(1,0.5,2),mua=rnorm(1),mus=rnorm(1),b0=rnorm(1),b=rnorm(20,0,1))}
parameters=c("b0","b0adj","b","badj","ba","baadj","bs","bsadj","tau","taua","taus","std","stda","stds")

cat("
model{
	for(i in 1:a){
		deaths[i]~dpois(lam[i])
		log(lam[i])<-log(years[i])+b0+b[i]+ba[age[i]]+bs[smoke[i]]
		b[i]~dnorm(0,tau)
		badj[i]<-b[i]-mean(b[])
	}
b0~dnorm(0,0.001)
b0adj<-b0+mean(b[])+mean(ba[])+mean(bs[])

for(ia in 1:5){
	ba[ia]~dnorm(mua,taua)
	baadj[ia]<-ba[ia]-mean(ba[])
}
mua~dnorm(0,0.001)

for(is in 1:4){
	bs[is]~dnorm(mus,taus)
	bsadj[is]<-bs[is]-mean(bs[])
}
mus~dnorm(0,0.001)

tau<-1/std/std
taua<-1/stda/stda
taus<-1/stds/stds
std~dunif(0,20)stda~dunif(0,20)stds~dunif(0,20)
}",file="death.txt")

deadsmokers.sim<-jags(data,inits,parameters,model.file="death.txt",n.chains=3,n.iter=20000,n.burnin=1000,n.thin=1)
attach(deadsmokers.sim)
traceplot(deadsmokers.sim)

write.csv(BUGSoutput$summary,"smokesummary.csv")

#b)
plot(deadsmokers.sim$BUGSoutput$sims.list$bs[56500:57000,4],type="l",main="bs[4]")
plot(deadsmokers.sim$BUGSoutput$sims.list$ba[56500:57000,4],type="l",main="ba[4]")
plot(deadsmokers.sim$BUGSoutput$sims.list$taua[56500:57000],type="l",main="taua")
autocorr.plot(deadsmokers.sim$BUGSoutput$sims.list$bs[,4],main="bs[4]",auto.layout=F)
autocorr.plot(deadsmokers.sim$BUGSoutput$sims.list$ba[,4],main="ba[4]",auto.layout=F)
autocorr.plot(deadsmokers.sim$BUGSoutput$sims.list$taua,main="taua",auto.layout=F)

#c)
prob <- exp(BUGSoutput$sims.list$bs[,4]-BUGSoutput$sims.list$bs[,1])
plot(density(prob))
mean(prob)
quantile(prob,c(0.025,0.975))
prob <- b[4]-b[1]

t(apply(deadsmokers.sim$BUGSoutput$sims.list$bs,2, 
    function(x){c(mean(x), sd(x),mean(x)/sd(x), quantile(x,probs=c(.025, .5,.975)))}  ))

#3.
#b)
    
x<-BUGSoutput$sims.array[,1,"b0"]
y<-BUGSoutput$sims.array[,2,"b0"]
z<-BUGSoutput$sims.array[,3,"b0"]

x<-BUGSoutput$sims.array[,1,"b0adj"]
y<-BUGSoutput$sims.array[,2,"b0adj"]
z<-BUGSoutput$sims.array[,3,"b0adj"]

x<-BUGSoutput$sims.array[,1,"taus"]
y<-BUGSoutput$sims.array[,2,"taus"]
z<-BUGSoutput$sims.array[,3,"taus"]

plot(density(z),col="blue")
lines(density(y),col="green")
lines(density(x),col="orange")

#c)

c(mean(x),mean(y),mean(z))

mclist<-as.mcmc.list(deadsmokers.sim$BUGSoutput)

batchSE(mclist)
effectiveSize(mclist)
summary(mclist)

#d)
geweke.diag(mclist)

g <- matrix(NA, nrow=nvar(mclist), ncol=2)
	for (v in 1:nvar(mclist)) {
		g[v,] <- gelman.diag(mclist[,v])$psrf
	}
h <- as.mcmc.list(g)
gelman.diag(as.mcmc(g))



#other
plot(x[18500:19000],type="l")

boxplot(deadsmokers.sim$BUGSoutput$sims.list$bs)

mc<-as.mcmc(deadsmokers.sim)
attach(mc)

detach(package:R2jags)

t<-as.mcmc(taus)
traceplot(t)
autocorr.plot(t)
plot(density(t))

b<-as.mcmc(b0)
traceplot(b)
autocorr.plot(b)

bb<-as.mcmc(b0adj)
traceplot(bb)
autocorr.plot(bb)

