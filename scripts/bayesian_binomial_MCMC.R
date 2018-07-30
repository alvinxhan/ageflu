require(rjags)
require(coda)

# clear environment completely
rm(list=ls())

# model specification
model_string <- "model {
xv ~ dbin(pv,nv)
pv ~ dbeta(1,1)
}"

# run MCMC
dataset=list(nv=108, xv=92) # inputs 
initialisation=list(pv=0.01)
jagmod=jags.model(textConnection(model_string),data=dataset,inits=initialisation, n.chains=1)
update(jagmod, n.iter=100000)
posterior <- coda.samples(jagmod, c("pv"),n.iter=1000000,thin=1000)

HPDinterval(posterior)
effectiveSize(posterior)
summary(posterior)
