args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=', args)][1])
BASE <- normalizePath(file.path(dirname(script), '..'), winslash='/', mustWork=TRUE)
setwd(BASE)
suppressPackageStartupMessages({library(INLA); library(spdep)})
dir.create('results',showWarnings=FALSE)
dir.create('results/model_objects',showWarnings=FALSE)
dir.create('environment',showWarnings=FALSE)
capture.output(sessionInfo(),file='environment/postreview_sessionInfo.txt')
Z <- read.csv('data/model_panel_numeric.csv',check.names=FALSE,fileEncoding='UTF-8-BOM')
G <- inla.read.graph('data/neighbors_228.graph')
stopifnot(nrow(Z)==3192, sum(Z$cases)==5634, length(unique(Z$region_id))==228)
stopifnot(all(Z$E_agestd>0), all(is.finite(as.matrix(Z[sprintf('V%02d',1:19)]))))
adj <- lapply(G$nbs,as.integer)
class(adj)<-'nb';attr(adj,'region.id')<-as.character(seq_along(adj))
lw<-nb2listw(adj,style='W',zero.policy=TRUE)
pc_bym2<-list(prec=list(prior='pc.prec',param=c(.5,.01)),phi=list(prior='pc',param=c(.5,2/3)))
pc_prec<-list(prec=list(prior='pc.prec',param=c(.5,.01)))
CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=.001)
parts<-c('V02','V04','V05','V14','V19')
LC<-read.csv('data/land_composition_percent.csv',fileEncoding='UTF-8-BOM')
stopifnot(length(LC$sigungu_cd)==nrow(Z),all(LC$sigungu_cd==Z$sigungu_cd),all(LC$year==Z$year))
audit<-read.csv('results/imputation_audit.csv',fileEncoding='UTF-8-BOM')
specs<-list(
 BASE=list(family='nbinomial',iid=TRUE,drop=character()),
 DROP50=list(family='nbinomial',iid=TRUE,drop=audit$variable_id[audit$filled_fraction>.5]),
 DROP20=list(family='nbinomial',iid=TRUE,drop=audit$variable_id[audit$filled_fraction>.2]),
 NB_NO_IID=list(family='nbinomial',iid=FALSE,drop=character()),
 POIS_IID=list(family='poisson',iid=TRUE,drop=character()),
 LAND_ALR=list(family='nbinomial',iid=TRUE,drop=character(),delta=.001),
 LAND_ALR_SMALL=list(family='nbinomial',iid=TRUE,drop=character(),delta=.0001))
all_coef<-list();all_fit<-list();all_hyper<-list();all_resid<-list();all_scales<-list()
for(nm in names(specs)) {
 s<-specs[[nm]];D<-Z;vars<-setdiff(sprintf('V%02d',1:19),s$drop)
 scales<-setNames(rep(NA_real_,length(parts)),parts)
 if(!is.null(s$delta)) {
  C<-as.matrix(LC[c(parts,'remainder')]);zero<-C==0; nz<-rowSums(zero)
  C<-C*(100-nz*s$delta)/100;C[zero]<-s$delta
  stopifnot(max(abs(rowSums(C)-100))<1e-8,all(C>0))
  for(v in parts) {
   a<-log(C[,v]/C[,'remainder']);scales[v]<-sd(a);D[[v]]<-(a-mean(a))/sd(a)
   all_scales[[paste(nm,v)]]<-data.frame(model=nm,variable_id=v,mean_logratio=mean(a),sd_logratio=sd(a),delta_pp=s$delta)
  }
 }
 rhs<-paste(vars,collapse=' + ')
 rhs<-paste(rhs,"f(region_id,model='bym2',graph=G,scale.model=TRUE,constr=TRUE,hyper=pc_bym2)",
            "f(year_id,model='rw1',scale.model=TRUE,hyper=pc_prec)",sep=' + ')
 if(s$iid)rhs<-paste(rhs,"f(id_st,model='iid',hyper=pc_prec)",sep=' + ')
 f<-as.formula(paste('cases ~',rhs))
 set.seed(1008);t0<-Sys.time()
 cat('\nSTART ',nm,' ',format(t0),'\n');flush.console()
 obj<-paste0('results/model_objects/',nm,'.rds')
 if(file.exists(obj)) {m<-readRDS(obj);cat('Using saved fit\n')} else {
  m<-inla(f,family=s$family,data=D,E=D$E_agestd,control.fixed=CF,
    control.compute=list(dic=TRUE,waic=TRUE,config=TRUE),control.predictor=list(compute=TRUE),num.threads='4:1')
  saveRDS(m,obj)
 }
 fx<-m$summary.fixed[vars,,drop=FALSE]
 cf<-data.frame(model=nm,variable_id=vars,beta=fx$mean,IRR=exp(fx$mean),lo=exp(fx$`0.025quant`),hi=exp(fx$`0.975quant`),
    credible=fx$`0.025quant`>0|fx$`0.975quant`<0,scale=ifelse(vars %in% parts & !is.null(s$delta),'standardized ALR','frozen selected form'),
    doubling_IRR=NA_real_,doubling_lo=NA_real_,doubling_hi=NA_real_)
 if(!is.null(s$delta))for(v in parts){j<-match(v,vars);mult<-log(2)/scales[v]
   cf[j,c('doubling_IRR','doubling_lo','doubling_hi')]<-exp(mult*unlist(fx[j,c('mean','0.025quant','0.975quant')]))}
 all_coef[[nm]]<-cf
 hp<-m$summary.hyperpar
 all_hyper[[nm]]<-data.frame(model=nm,parameter=rownames(hp),mean=hp$mean,sd=hp$sd,lo=hp$`0.025quant`,median=hp$`0.5quant`,hi=hp$`0.975quant`)
 size<-if(s$family=='nbinomial')hp$mean[grep('size',rownames(hp),ignore.case=TRUE)][1] else Inf

 mu<-D$E_agestd*m$summary.fitted.values$mean[seq_len(nrow(D))]
 stopifnot(all(is.finite(mu)),sum(mu)/sum(D$cases)>.5,sum(mu)/sum(D$cases)<1.5)
 pr<-tapply((D$cases-mu)/sqrt(mu+mu^2/size),D$region_id,mean)
 set.seed(1008);mc<-moran.mc(as.numeric(pr),lw,nsim=999,zero.policy=TRUE)
 all_resid[[nm]]<-data.frame(model=nm,region_id=seq_along(pr),pearson_residual=as.numeric(pr))
 all_fit[[nm]]<-data.frame(model=nm,family=s$family,n=nrow(D),p=length(vars),DIC=m$dic$dic,WAIC=m$waic$waic,
   nb_size=size,pearson_moran=as.numeric(mc$statistic),pearson_p=mc$p.value,fitted_cases=sum(mu),
   elapsed_seconds=as.numeric(difftime(Sys.time(),t0,units='secs')))
 write.csv(do.call(rbind,all_coef),'results/sensitivity_coefficients.csv',row.names=FALSE)
 write.csv(do.call(rbind,all_fit),'results/sensitivity_fit.csv',row.names=FALSE)
 write.csv(do.call(rbind,all_hyper),'results/hyperparameters.csv',row.names=FALSE)
 write.csv(do.call(rbind,all_resid),'results/pearson_residuals.csv',row.names=FALSE)
 if(length(all_scales))write.csv(do.call(rbind,all_scales),'results/alr_scaling.csv',row.names=FALSE)
 cat(sprintf('DONE %s WAIC=%.2f Moran=%.3f p=%.3f (%0.fs)\n',nm,m$waic$waic,mc$statistic,mc$p.value,all_fit[[nm]]$elapsed_seconds));flush.console()
}
cat('\nCOMPLETE: all seven locked specifications fitted.\n')
