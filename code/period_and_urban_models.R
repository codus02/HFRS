args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=', args)][1])
BASE <- normalizePath(file.path(dirname(script), '..'),winslash='/',mustWork=TRUE)
setwd(BASE)
suppressPackageStartupMessages(library(INLA))
dir.create('results/secondary',showWarnings=FALSE,recursive=TRUE)
Z <- read.csv('data/model_panel_numeric.csv',check.names=FALSE)
U <- read.csv('data/urban_sensitivity_covariate.csv',check.names=FALSE)
stopifnot(all(U$sigungu_cd==Z$sigungu_cd),all(U$year==Z$year))
Z$U01 <- U$U01_urban_population_share_standardized
G <- inla.read.graph('data/neighbors_228.graph')
pc_bym2<-list(prec=list(prior='pc.prec',param=c(.5,.01)),phi=list(prior='pc',param=c(.5,2/3)))
pc_prec<-list(prec=list(prior='pc.prec',param=c(.5,.01)))
CF<-list(mean=0,prec=4,mean.intercept=0,prec.intercept=.001)
sp <- "f(region_id,model='bym2',graph=G,scale.model=TRUE,constr=TRUE,hyper=pc_bym2)"
ic <- "f(region_id,model='besag',graph=G,scale.model=TRUE,constr=TRUE,hyper=pc_prec)"
ti <- "f(year_id,model='rw1',scale.model=TRUE,hyper=pc_prec)"
it <- "f(id_st,model='iid',hyper=pc_prec)"
specs<-list(P1=list(years=2011:2015,extra=character(),re=c(sp,ti,it)),
            P2=list(years=2016:2019,extra=character(),re=c(sp,ti,it)),
            P3=list(years=2020:2024,extra=character(),re=c(sp,ti,it)),
            URBAN=list(years=2011:2024,extra='U01',re=c(sp,ti,it)))
if('--ladder' %in% commandArgs(trailingOnly=TRUE)) {
 re<-list(M1=character(),M2=ic,M3=sp,M4=c(sp,it),M5=c(sp,ti),M6=c(sp,ti,it))
 for(n in names(re))specs[[n]]<-list(years=2011:2024,extra=character(),re=re[[n]])
}
coef<-list();fits<-list();spatial<-list()
for(nm in names(specs)) {
 s<-specs[[nm]];D<-Z[Z$year %in% s$years,];D<-D[order(D$region_id,D$year),]
 D$year_id<-match(D$year,sort(unique(D$year)));D$id_st<-seq_len(nrow(D))
 vars<-c(sprintf('V%02d',1:19),s$extra)
 f<-as.formula(paste('cases ~',paste(c(vars,s$re),collapse=' + ')))
 set.seed(1008);cat('START ',nm,'\n');flush.console()
 m<-inla(f,family='nbinomial',data=D,E=D$E_agestd,control.fixed=CF,
          control.compute=list(dic=TRUE,waic=TRUE),control.predictor=list(compute=TRUE),num.threads='4:1')
 b<-m$summary.fixed[vars,]
 coef[[nm]]<-data.frame(model=nm,variable_id=vars,IRR=exp(b$mean),lo=exp(b$`0.025quant`),hi=exp(b$`0.975quant`))
 fits[[nm]]<-data.frame(model=nm,n=nrow(D),cases=sum(D$cases),DIC=m$dic$dic,WAIC=m$waic$waic)
 if(any(s$re==sp)) {
  sr<-m$summary.random$region_id;ids<-unique(D[c('region_id','sigungu_cd')])

  sr<-sr[sr$ID<=G$n,]
  spatial[[nm]]<-data.frame(model=nm,region_id=sr$ID,sigungu_cd=ids$sigungu_cd[match(sr$ID,ids$region_id)],
                            mean=sr$mean,lo=sr$`0.025quant`,hi=sr$`0.975quant`,spatial_multiplier=exp(sr$mean))
 }
 write.csv(do.call(rbind,coef),'results/secondary/coefficients.csv',row.names=FALSE)
 write.csv(do.call(rbind,fits),'results/secondary/fit.csv',row.names=FALSE)
 write.csv(do.call(rbind,spatial),'results/secondary/spatial.csv',row.names=FALSE)
 cat('DONE ',nm,' WAIC ',m$waic$waic,'\n');flush.console()
}
