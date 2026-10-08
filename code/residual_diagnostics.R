args<-commandArgs(trailingOnly=FALSE)
script<-sub('^--file=','',args[grepl('^--file=',args)][1])
BASE<-normalizePath(file.path(dirname(script),'..'),winslash='/');setwd(BASE)
suppressPackageStartupMessages({library(INLA);library(spdep)})
Z<-read.csv('data/model_panel_numeric.csv');g<-inla.read.graph('data/neighbors_228.graph')
nb<-lapply(g$nbs,as.integer);class(nb)<-'nb';attr(nb,'region.id')<-as.character(1:228)
lw<-nb2listw(nb,style='W',zero.policy=TRUE)
m<-readRDS('results/model_objects/BASE.rds')
O<-as.numeric(tapply(Z$cases,Z$region_id,sum));E<-as.numeric(tapply(Z$E_agestd,Z$region_id,sum))
set.seed(1008);raw<-moran.mc(O/E,lw,nsim=999,zero.policy=TRUE)

mu<-m$summary.fitted.values$mean*Z$E_agestd
sz<-m$summary.hyperpar[grep('size',rownames(m$summary.hyperpar)), 'mean']
pr<-as.numeric(tapply((Z$cases-mu)/sqrt(mu+mu^2/sz),Z$region_id,mean))
set.seed(1008);pearson<-moran.mc(pr,lw,nsim=999,zero.policy=TRUE)
set.seed(1008)
sm<-inla.posterior.sample(1000,m,seed=1008,num.threads='1:1')
idx<-grep('^Predictor:',rownames(sm[[1]]$latent));stopifnot(length(idx)==nrow(Z))
ids<-as.integer(sub('^Predictor:','',rownames(sm[[1]]$latent)[idx]));idx<-idx[order(ids)]
rep<-vapply(sm,function(s){
 eta<-s$latent[idx,1];size<-s$hyperpar[grep('size',names(s$hyperpar),ignore.case=TRUE)]
 as.numeric(tapply(rnbinom(nrow(Z),size=size,mu=exp(eta)*Z$E_agestd),Z$region_id,sum))
},numeric(228))
set.seed(1008)
u<-vapply(1:228,function(i)runif(1,mean(rep[i,]<O[i]),mean(rep[i,]<=O[i])),numeric(1))
qr<-qnorm(pmin(pmax(u,.0005),.9995))
set.seed(1008);quantile<-moran.mc(qr,lw,nsim=999,zero.policy=TRUE)
out<-data.frame(diagnostic=c('Observed SIR','Posterior predictive quantile','Pearson (offset corrected)'),
 Moran_I=c(raw$statistic,quantile$statistic,pearson$statistic),p=c(raw$p.value,quantile$p.value,pearson$p.value),
 permutations=999,posterior_draws=c(NA,1000,NA))
write.csv(out,'results/corrected_residual_diagnostics.csv',row.names=FALSE)
write.csv(data.frame(region_id=1:228,observed=O,expected=E,quantile_residual=qr,pearson_residual=pr),
          'results/corrected_district_residuals.csv',row.names=FALSE)

sdrows<-lapply(grep('Precision',names(m$marginals.hyperpar),value=TRUE),function(nm){
 marg<-inla.tmarginal(function(tau)1/sqrt(tau),m$marginals.hyperpar[[nm]])
 data.frame(component=nm,posterior_mean_sd=inla.emarginal(identity,marg),
 median=inla.qmarginal(.5,marg),lo=inla.qmarginal(.025,marg),hi=inla.qmarginal(.975,marg))
})
write.csv(do.call(rbind,sdrows),'results/random_effect_sd.csv',row.names=FALSE)
print(out);cat('Corrected diagnostics completed.\n')
