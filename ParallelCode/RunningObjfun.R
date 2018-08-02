library(glmm)
library(parallel)
library(doParallel)
library(foreach)
library(itertools)

set.seed(1234)

data("BoothHobert")
test<-glmm(y~0+x1,list(y~0+z1),varcomps.names=c("z1"),data=BoothHobert,
          family.glmm=bernoulli.glmm,m=100,doPQL=FALSE,debug=TRUE)

mod.mcml<-test$mod.mcml
debug<-test$debug
nu.pql<-debug$nu.pql
beta.pql<-debug$beta.pql
family.glmm<-test$family.glmm
umat<-debug$umat
u.pql<-debug$u.star
m1<-debug$m1
ntrials<-1

par<-c(6,1.5)
del<-rep(10^-8,2)

objfun<-glmm:::objfun
getEk<-glmm:::getEk
addVecs<-glmm:::addVecs

eek<-getEk(mod.mcml$z)
Aks<-Map("*",eek,nu.pql)
D.star<-addVecs(Aks) 
D.star<-diag(D.star)
D.star.inv<-solve(D.star)

Z=do.call(cbind,mod.mcml$z)
eta.star<-as.vector(mod.mcml$x%*%beta.pql+Z%*%u.pql)
cdouble<-bernoulli.glmm()$cpp(eta.star) #still a vector
cdouble<-diag(cdouble)
Sigmuh.inv<- t(Z)%*%cdouble%*%Z+D.star.inv
Sigmuh<-solve(Sigmuh.inv)

p1=p2=p3=1/3
zeta=5

tconstant<-glmm:::tconstant
no_cores <- max(1, detectCores() - 1)
getFamily<-glmm:::getFamily

regobj<-objfun(par=par, nbeta=1, nu.pql=nu.pql, umat=umat, u.star=u.pql, mod.mcml=mod.mcml, family.glmm=family.glmm,p1=p1,p2=p2,p3=p3,m1=m1, Sigmuh=Sigmuh, D.star=D.star, Sigmuh.inv= Sigmuh.inv, zeta=zeta, ntrials=ntrials)

objfun2 <-
  function(par, nbeta, nu.pql, umat, u.star, mod.mcml, family.glmm, cache, p1, p2, p3, m1, D.star, Sigmuh, Sigmuh.inv, zeta, ntrials){
    
    beta<-par[1:nbeta]
    nu<-par[-(1:nbeta)]
    m<-nrow(umat)
    
    if (!missing(cache)) stopifnot(is.environment(cache))
    
    if(any(nu<=0)){
      out<-list(value=-Inf,gradient=rep(1,length(par)),hessian=as.matrix(c(rep(1,length(par)^2)),nrow=length(par)))
      return(out)
    }
    
    Z=do.call(cbind,mod.mcml$z)
    T<-length(mod.mcml$z)
    nrand<-lapply(mod.mcml$z,ncol)
    nrandom<-unlist(nrand)
    
    
    
    family.glmm<-getFamily(family.glmm)
    if(family.glmm$family.glmm=="bernoulli.glmm"){family_glmm=1}	
    if(family.glmm$family.glmm=="poisson.glmm"){family_glmm=2}	
    if(family.glmm$family.glmm=="binomial.glmm"){family_glmm=3}	
    
    Dstarinvdiag<-1/diag(D.star)
    D.star.inv<-diag(Dstarinvdiag)
    
    logdet.D.star.inv<-	-sum(log(diag(D.star)))
    logdet.Sigmuh.inv<-sum(log(eigen(Sigmuh.inv,symmetric=TRUE)$values))
    myq<-nrow(D.star.inv)
    
    tconst<-tconstant(zeta,myq,Dstarinvdiag)
    
    #for the particular value of nu we're interested in, need to prep for distRandGenC
    eek<-getEk(mod.mcml$z)
    preDinvfornu<-Map("*",eek,(1/nu))
    Dinvfornu<-addVecs(preDinvfornu)
    logdetDinvfornu<-sum(log(Dinvfornu))
    Dinvfornu <- diag(Dinvfornu)
    
    meow<-rep(1,T+1)
    meow[1]<-0
    throwaway<-T+1
    meow[2:throwaway]<-cumsum(nrandom)
    
    pea<-c(p1,p2,p3)
    n <- nrow(mod.mcml$x)
    
    ##need to scale first m1 vectors of generated random effects by multiplying by A
    
    #	preAfornu<-Map("*",eek,sqrt(nu))
    #	Afornu<-addVecs(preAfornu)
    
    #	for(k in 1:m1){
    #		u.swoop<-umat[k,]
    #		umat[k,]<-u.swoop*Afornu
    #		}
    
    stuff<-.C("objfunc", as.double(mod.mcml$y),as.double(t(umat)), as.integer(myq), as.integer(m), as.double(mod.mcml$x), as.integer(n), as.integer(nbeta), as.double(beta), as.double(Z), as.double(Dinvfornu), as.double(logdetDinvfornu),as.integer(family_glmm), as.double(D.star.inv), as.double(logdet.D.star.inv), as.double(u.star), as.double(Sigmuh.inv), as.double(logdet.Sigmuh.inv), pea=as.double(pea), nps=as.integer(length(pea)), T=as.integer(T), nrandom=as.integer(nrandom), meow=as.integer(meow),nu=as.double(nu), zeta=as.integer(zeta),tconst=as.double(tconst), v=double(m), ntrials=as.integer(ntrials), value=double(1),gradient=double(length(par)),hessian=double((length(par))^2))
    
    #R CMD BUILD name
    #R CMD check tar --as-cran
    #--> description change version 1.2.4, change date, add contributer
    
    
    if (!missing(cache)) cache$weights<-stuff$v		
    
    list(stuff)
    
    list(value=stuff$value,gradient=stuff$gradient,hessian=matrix(stuff$hessian,ncol=length(par),byrow=FALSE))
    
  }

newobj<-objfun(par=par, nbeta=1, nu.pql=nu.pql, umat=umat, u.star=u.pql, mod.mcml=mod.mcml, family.glmm=family.glmm,p1=p1,p2=p2,p3=p3,m1=m1, Sigmuh=Sigmuh, D.star=D.star, Sigmuh.inv= Sigmuh.inv, zeta=zeta, ntrials=ntrials, no_cores=no_cores)

sum(newobj[[4]][[1]]$hessian[1], newobj[[4]][[2]]$hessian[1], newobj[[4]][[3]]$hessian[1])

newobj[[1]][[1]]$value
newobj[[1]][[2]]$value
newobj[[1]][[3]]$value

mean(c(newobj[[1]][[1]]$value, newobj[[1]][[2]]$value, newobj[[1]][[3]]$value))

newobj[[1]][[1]]$gradient[1]
newobj[[1]][[2]]$gradient
newobj[[1]][[3]]$gradient

mean(c(newobj[[1]][[1]]$gradient[1], newobj[[1]][[2]]$gradient[1], newobj[[1]][[3]]$gradient[1]))

newobj[[1]][[1]]$gradient[1]+ newobj[[1]][[2]]$gradient[1]+ newobj[[1]][[3]]$gradient[1]



res <- c()
for(i in 1:no_cores){
  res[i] <- newobj[[1]][[i]][[4]]*exp(newobj[[1]][[i]]$value)
  expadd <- expadd + res[i]
}

hessmat <- matrix(rep(0, 12), nrow = 4, ncol = 3)
hess <- c(rep(0, 4))
for(j in 1:length(newobj[[1]][[1]]$hessian)){
  hessadd <- 0
  for(i in 1:length(res)){
    hessmat[j,i] <- res[i]*newobj[[1]][[i]]$hessian[j]
    hessadd <- hessadd + res[i]*newobj[[1]][[i]]$hessian[j]
  }
  hess[j] <- hessadd/expadda
}

grad <- matrix(rep(0, 6), nrow = 2, ncol = 3)
for(i in 1:length(res)){
  gradients[[i]] <- (res[i]/expadd)*gradients[[i]]
}
sum(grad[1,])

val <- 0
for(i in 1:no_cores){
  val <- val + sum(newobj[[1]][[i]]$v)
}
log(val/m)

vals <- c()
for(i in 1:no_cores){
  vals[i] <- newobj[[1]][[i]]$value
}
a <- max(vals)

res <- c(rep(0, no_cores))
expadda <- 0
for(i in 1:no_cores){
  res[i] <- newobj[[1]][[i]][[4]]*exp(newobj[[1]][[i]]$value - a)
  expadda <- expadda + newobj[[1]][[i]][[4]]*exp(newobj[[1]][[i]]$value - a)
}

hessian <- c(rep(0, length(newobj[[4]][[1]]$hessian)))
for(j in 1:length(newobj[[4]][[1]]$hessian)){
  hessadd <- 0
  for(i in 1:no_cores){
    hessadd <- hessadd + newobj[[4]][[i]]$hessian[j]
  }
  hessian[j] <- hessadd
}