args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=', args)][1])
BASE <- normalizePath(file.path(dirname(script), '..'), winslash='/', mustWork=TRUE)
setwd(BASE)
suppressPackageStartupMessages({library(INLA); library(spdep); library(sf)})
dir.create('results/extended', showWarnings=FALSE, recursive=TRUE)
dir.create('results/model_objects', showWarnings=FALSE)
Z <- read.csv('data/model_panel_numeric.csv', check.names=FALSE, fileEncoding='UTF-8-BOM')
G0 <- inla.read.graph('data/neighbors_228.graph')
stopifnot(nrow(Z)==3192, sum(Z$cases)==5634, length(unique(Z$region_id))==228)
vars_all <- sprintf('V%02d', 1:19)
key <- unique(Z[c('sigungu_cd','region_id')]); key <- key[order(key$region_id),]
stopifnot(all(key$region_id==1:228))

geo <- st_read('data/district_boundaries.gpkg', quiet=TRUE)
geo <- geo[match(key$sigungu_cd, geo$sigungu_cd),]; stopifnot(nrow(geo)==228, !anyNA(geo$sigungu_cd))
geo <- st_make_valid(geo)
cent <- st_coordinates(st_centroid(st_geometry(geo)))
supplied_nb <- lapply(G0$nbs, as.integer); class(supplied_nb) <- 'nb'; attr(supplied_nb,'region.id') <- as.character(1:228)
connect_components <- function(nb, note) {

  added <- character()
  for (i in which(card(nb)==0)) { nb[[i]] <- sort(unique(c(nb[[i]][nb[[i]]>0], supplied_nb[[i]]))); for (j in supplied_nb[[i]]) nb[[j]] <- sort(unique(c(nb[[j]][nb[[j]]>0], i))); added <- c(added, sprintf('island %d <- supplied links %s', i, paste(supplied_nb[[i]], collapse='/'))) }
  repeat {
    comp <- n.comp.nb(nb); if (comp$nc==1) break
    main <- which.max(tabulate(comp$comp.id)); small <- which(comp$comp.id!=main); big <- which(comp$comp.id==main)
    d <- as.matrix(dist(cent))[small, big, drop=FALSE]; w <- which(d==min(d), arr.ind=TRUE)[1,]
    a <- small[w[1]]; b <- big[w[2]]
    nb[[a]] <- sort(unique(c(nb[[a]][nb[[a]]>0], b))); nb[[b]] <- sort(unique(c(nb[[b]][nb[[b]]>0], a)))
    added <- c(added, sprintf('component link %d-%d (%.1f km)', a, b, min(d)/1000))
  }
  list(nb=nb, log=paste(c(note, added), collapse='; '))
}
graphs <- list(SUPPLIED=list(nb=supplied_nb, log='supplied graph (592 undirected edges)'))
q <- poly2nb(geo, queen=TRUE); graphs$QUEEN <- connect_components(q, sprintf('queen contiguity, %d islands', sum(card(q)==0)))
r <- poly2nb(geo, queen=FALSE); graphs$ROOK <- connect_components(r, sprintf('rook contiguity, %d islands', sum(card(r)==0)))
for (k in 2:7) { nbk <- make.sym.nb(knn2nb(knearneigh(cent, k=k))); graphs[[paste0('KNN', k)]] <- connect_components(nbk, sprintf('%d nearest neighbours, symmetrised', k)) }
graph_summary <- do.call(rbind, lapply(names(graphs), function(g) { nb <- graphs[[g]]$nb; data.frame(graph=g, edges=sum(card(nb))/2, mean_neighbours=mean(card(nb)), components=n.comp.nb(nb)$nc, construction=graphs[[g]]$log) }))
write.csv(graph_summary, 'results/extended/graph_summary.csv', row.names=FALSE)
to_graph <- function(nb) { f <- tempfile(fileext='.graph'); nb2INLA(f, nb); inla.read.graph(f) }

CF0 <- list(mean=0, prec=4, mean.intercept=0, prec.intercept=.001)
prior_set <- function(u=.5, a=.01, phi_u=.5, phi_a=2/3) list(bym2=list(prec=list(prior='pc.prec', param=c(u,a)), phi=list(prior='pc', param=c(phi_u,phi_a))), prec=list(prec=list(prior='pc.prec', param=c(u,a))))
fit_one <- function(nm, D, graph, family='nbinomial', vars=vars_all, priors=prior_set(), CF=CF0, iid=TRUE, cpo=FALSE, extra='') {
  rhs <- paste(c(vars, extra[nzchar(extra)]), collapse=' + ')
  rhs <- paste(rhs, "f(region_id,model='bym2',graph=graph,scale.model=TRUE,constr=TRUE,hyper=priors$bym2)",
               "f(year_id,model='rw1',scale.model=TRUE,hyper=priors$prec)", sep=' + ')
  if (iid) rhs <- paste(rhs, "f(id_st,model='iid',hyper=priors$prec)", sep=' + ')
  f <- as.formula(paste('cases ~', rhs))
  obj <- paste0('results/model_objects/ext_', nm, '.rds'); t0 <- Sys.time()
  cat('\nSTART', nm, format(t0), '\n'); flush.console()
  if (file.exists(obj)) { m <- readRDS(obj); cat('Using saved fit\n') } else {
    set.seed(1008)
    m <- inla(f, family=family, data=D, E=D$E_agestd, control.fixed=CF,
              control.compute=list(dic=TRUE, waic=TRUE, cpo=cpo, config=FALSE), control.predictor=list(compute=TRUE), num.threads='4:1')
    saveRDS(m, obj)
  }
  fx <- m$summary.fixed[c(vars, extra[nzchar(extra)]),,drop=FALSE]
  cf <- data.frame(model=nm, variable_id=rownames(fx), beta=fx$mean, IRR=exp(fx$mean), lo=exp(fx$`0.025quant`), hi=exp(fx$`0.975quant`), credible=fx$`0.025quant`>0 | fx$`0.975quant`<0)
  hp <- m$summary.hyperpar
  hy <- data.frame(model=nm, parameter=rownames(hp), mean=hp$mean, sd=hp$sd, lo=hp$`0.025quant`, hi=hp$`0.975quant`)
  mu <- D$E_agestd * m$summary.fitted.values$mean[seq_len(nrow(D))]
  ft <- data.frame(model=nm, family=family, n=nrow(D), p=length(vars)+sum(nzchar(extra)), DIC=m$dic$dic, WAIC=m$waic$waic, fitted_cases=sum(mu), observed_cases=sum(D$cases),
                   elapsed_seconds=as.numeric(difftime(Sys.time(), t0, units='secs')))
  cat(sprintf('DONE %s WAIC=%.2f credible=%d (%0.fs)\n', nm, m$waic$waic, sum(cf$credible), ft$elapsed_seconds)); flush.console()
  list(coef=cf, hyper=hy, fit=ft, model=m, mu=mu)
}
COEF <- list(); HYP <- list(); FIT <- list()
record <- function(res) { COEF[[res$fit$model]] <<- res$coef; HYP[[res$fit$model]] <<- res$hyper; FIT[[res$fit$model]] <<- res$fit
  write.csv(do.call(rbind, COEF), 'results/extended/coefficients.csv', row.names=FALSE)
  write.csv(do.call(rbind, HYP), 'results/extended/hyperparameters.csv', row.names=FALSE)
  write.csv(do.call(rbind, FIT), 'results/extended/fit.csv', row.names=FALSE) }

res <- fit_one('BASE_CPO', Z, G0, cpo=TRUE); record(res)
m <- res$model
pit <- data.frame(sigungu_cd=Z$sigungu_cd, year=Z$year, cases=Z$cases, expected=Z$E_agestd, fitted_cases=res$mu, pit=m$cpo$pit, cpo=m$cpo$cpo, cpo_failure=m$cpo$failure)
write.csv(pit, 'results/extended/pit_cpo.csv', row.names=FALSE)
agg <- aggregate(cbind(cases, E_agestd, fitted=res$mu) ~ sigungu_cd + region_id, data=cbind(Z, fitted=res$mu), FUN=sum)
agg$observed_SIR <- agg$cases/agg$E_agestd; agg$predicted_SIR <- agg$fitted/agg$E_agestd; agg$log_ratio_obs_pred <- log((agg$cases+.5)/(agg$fitted+.5))
write.csv(agg, 'results/extended/district_observed_predicted.csv', row.names=FALSE)

for (g in setdiff(names(graphs), 'SUPPLIED')) record(fit_one(paste0('GRAPH_', g), Z, to_graph(graphs[[g]]$nb)))

record(fit_one('PRIOR_PCPREC_1_001', Z, G0, priors=prior_set(1, .01)))
record(fit_one('PRIOR_PCPREC_1_005', Z, G0, priors=prior_set(1, .05)))
record(fit_one('PRIOR_PCPREC_3_005', Z, G0, priors=prior_set(3, .05)))
record(fit_one('PRIOR_PHI_050', Z, G0, priors=prior_set(.5, .01, .5, .5)))
record(fit_one('PRIOR_PHI_033', Z, G0, priors=prior_set(.5, .01, .5, 1/3)))
record(fit_one('PRIOR_FIXED_DIFFUSE', Z, G0, CF=list(mean=0, prec=.001, mean.intercept=0, prec.intercept=.001)))

record(fit_one('ZINB_IID', Z, G0, family='zeroinflatednbinomial1'))

Zx <- Z[!(Z$year %in% 2020:2022),]; record(fit_one('EXCL_2020_2022', Zx, G0))
Zi <- Z; Zi$covid_2020_2022 <- as.integer(Zi$year %in% 2020:2022); record(fit_one('COVID_INDICATOR', Zi, G0, extra='covid_2020_2022'))

record(fit_one('DROP_HEALTHCARE', Z, G0, vars=setdiff(vars_all, c('V07','V09','V18'))))

X <- as.matrix(Z[vars_all])
vif <- sapply(vars_all, function(v) { r2 <- summary(lm(X[,v] ~ X[,setdiff(vars_all, v)]))$r.squared; 1/(1-r2) })
rho <- cor(X, method='spearman'); pairs <- which(abs(rho)>.7 & upper.tri(rho), arr.ind=TRUE)
coll <- data.frame(diagnostic=c('Covariates, no.','Maximum variance inflation factor','Covariates with VIF > 5, no.','Covariates with VIF > 10, no.','Design-matrix condition number (scaled)','Covariate pairs with |Spearman rho| > 0.7, no.'),
                   value=c(19, round(max(vif),2), sum(vif>5), sum(vif>10), round(kappa(scale(X), exact=TRUE),1), nrow(pairs)))
write.csv(coll, 'results/extended/collinearity_summary.csv', row.names=FALSE)
write.csv(data.frame(variable_id=vars_all, vif=round(vif,3)), 'results/extended/vif_by_covariate.csv', row.names=FALSE)
if (nrow(pairs)) write.csv(data.frame(a=vars_all[pairs[,1]], b=vars_all[pairs[,2]], spearman=round(rho[pairs],3)), 'results/extended/high_correlation_pairs.csv', row.names=FALSE)

base <- read.csv('results/sensitivity_coefficients.csv'); base <- base[base$model=='BASE',]
cred <- base$variable_id[base$credible]
cf <- do.call(rbind, COEF); gm <- cf[grepl('^GRAPH_', cf$model),]
rob <- do.call(rbind, lapply(cred, function(v) { d <- gm[gm$variable_id==v,]; d$graph <- sub('GRAPH_','',d$model)
  n <- sum(d$credible); contig <- all(d$credible[d$graph %in% c('QUEEN','ROOK')]); knn <- sum(d$credible[grepl('KNN', d$graph)])
  data.frame(variable_id=v, reference_IRR=base$IRR[base$variable_id==v], irr_min=min(d$IRR), irr_max=max(d$IRR), graphs_credible=n,
             queen=d$credible[d$graph=='QUEEN'], rook=d$credible[d$graph=='ROOK'], knn_credible_of_6=knn,
             class=ifelse(n==8, 'Graph-robust', ifelse(contig & knn<=2, 'Contiguity-only', ifelse(n>=5, 'Mostly robust', 'Graph-sensitive')))) }))
write.csv(rob, 'results/extended/graph_robustness.csv', row.names=FALSE)
capture.output(sessionInfo(), file='environment/extended_sessionInfo.txt')
cat('\nCOMPLETE: extended sensitivity program finished.\n')
