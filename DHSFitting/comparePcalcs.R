rm(list=ls())

## jobnum=1;group.ind=1;s.epic=1;s.demog=1;s.bmb=1;s.bfb=1;s.bme=1;s.bfe=1;s.bmp=1;s.bfp=1;acute.sc=1;bmp.sc=1;bfp.sc=1;late.sc=1;aids.sc=1;new.folder=TRUE;simul=FALSE;maxN=1000;seed.bump=1;sigma.ad=TRUE;batchdirnm="results/DHSFits/BurundiRealAc";psNonPar=TRUE;infl.fac=50;death=TRUE;bmb.sc=1;bfb.sc=1;bme.sc=1;bfe.sc=1;het.b=FALSE;het.b.sd=0;het.b.cor=0;het.e=FALSE;het.e.sd=0;het.e.cor=0;het.p=FALSE;het.p.sd=0;het.p.cor=0;het.gen=FALSE;het.gen.sd=0;het.gen.cor=0;het.beh=FALSE;het.beh.sd=0;het.beh.cor=0;scale.by.sd=TRUE;scale.adj=1;sample.tmar=FALSE;tmar=(65*12):(111*12);each=200;tint=111*12;short.test=F;all.cores=T;nc=6;adapt=TRUE;d.nburn=400;d.nthin=3;d.niter=1200;nburn=500;nthin=1;niter=8000;survive=T;tell=100;low.coverage.arv=F;partner.arv=F;fsd.sens=F

## setwd('/home1/02413/sbellan/DHSProject/DHSFitting/')
## ## This script must be called by R CMD BATCH with a set of arguments. See DHSFitMK.R which makes
## ## DHSFitControlFile.txt, the latter has one line with a set of arguments to be sent to a cluster.
## args <- commandArgs(TRUE)
## ## args is now a list of character vectors, we cycle through each element of the list and evaluate the expressions.
## if(length(args)>0) { for(ii in 1:length(args))  eval(parse(text=args[[ii]])) }
## source('../SDPSimulations/SimulationFunctions.R') # simulating functions (from other project folder)
## load("data files/copula sigmas.Rdata")

## ##load("data files/alldhs.Rdata")         # DHS data
## load("data files/allDHSAIS.Rdata")         # DHS data
## load("data files/ds.nm.all.Rdata") # country names
## load("data files/pars.arr.ac.Rdata")       # previously fitted parameters
## load("data files/csurv.Rdata")             # AIDS survival times
## ##  [1] "Congo"      "DRC"        "Ethiopia"   "Kenya"      "Lesotho"   
## ##  [6] "Malawi"     "Mozambique" "Rwanda"     "Swaziland"  "Tanzania"  
## ## [11] "Uganda"     "WA"         "Zambia"     "Zimbabwe"  
## if(low.coverage.arv) { # if we're assuming 50% of ARV coverage results in no transmission
##     load("data files/allepicm.5.Rdata") # deprecated (need to reload these files)
##     load("data files/allepicf.5.Rdata")
##   }else{ # if we're assuming 100% of ARV coverage results in no transmission
##     load("data files/epic.Rdata")       # loads epicm & epicf
##   }
## hazs <- c("bmb","bfb","bme","bfe","bmp","bfp") ## transmission coefficient names
## parnames <- c("bmb","bfb","bme","bfe","bmp","lrho") # parameters to fit
## trans.ratio <- 1 ## m->f : f<-m; within couples
## library(ade4); library(mvtnorm);library(mnormt);library(multicore);library(coda);library(abind); library(plotrix);
## ####################################################################################################
## ## West Africa is pooled and analyzed in a separate script modified to
## ## deal with multiple countrie's prevalences simultaneously.
## ######################################################################
## print(batchdirnm) ## where are we saving things?

## ## These parameters will be used for simulating data that are then fit.
## simpars <- c(bmb = .01, bfb = .03, bme = .012, bfe =.01, bmp=.03, lrho = 0)
## simpars <- c(simpars, bfp = as.numeric(simpars['bmp'] * exp(simpars['lrho']))) # bfp as a function of bmp & rho
## ## Set up new directory to store results
## group <- levels(dat$group)[group.ind]
## dirnm <- paste0(group,'RealAc',acute.sc, '-jb', jobnum)
## step <- 1
## ndirnm <- paste0(dirnm,'-', step)
## while(new.folder & file.exists(ndirnm)) { ## don't continue from old analysis (otherwise redo new folder)
##     ndirnm <- paste0(dirnm,'-', step)
##     step <- step + 1
## }
## ndirnm <- file.path(batchdirnm, ndirnm)
## if(!file.exists(ndirnm))        dir.create(ndirnm) ## create directory
## wa <- levels(dat$group)[group.ind]=="WA"
## if(!wa) { ## if not WA just set countries
##     group <- levels(dat$group)[group.ind]
##     print(paste(nlevels(dat$ds),"country data sets"))
##     ds <- levels(dat$ds)[grepl(group, levels(dat$ds))]   # data set we're working on
##     dat <- dat[dat$ds %in% ds,]           # only looking at that data set
##     print(paste("analyzing",paste(ds, collapse = " & "),"couples data using", dat$epic.nm[1],
##                 "epidemic curve"))
##     if(fsd.sens) { # if doing sensitivity analysis to female sexual debut
##         ## Of females saying their sexual debut occured at
##         ## marriage, assume 30% were lying and that actually
##         ## sexual debut occurred one year earlier.
##         print(paste("lowering female sexual debut by one year for 30% of females that stated they first started having sex at marriage (sensitivity analysis)"))
##         fsd.mar.ind <- which(dat$tfs==dat$tmar) #which sexual debuts occurred at marriage
##         alter.ind <- sample(fsd.mar.ind, size = round(length(fsd.mar.ind)*.3)) # ones to alter
##         dat$tfs[alter.ind] <- dat$tfs[alter.ind] - 12
##     }
##     save(dat, file = file.path(ndirnm, paste0(group,".Rdata"))) # just to indicate which ds working on when looking in file dir
## }else{ ## For West Africa, need to select all countries in region
##     group <- levels(dat$group)[group.ind]
##     print(paste(nlevels(dat$ds),"country data sets"))
##     dat <- dat[dat$group == group,]                 # west africa
##     if(fsd.sens) { # if doing sensitivity analysis to female sexual debut
##         ## Of females saying their sexual debut occured at
##         ## marriage, assume 30% were lying and that actually
##         ## sexual debut occurred one year earlier.
##         print(paste("lowering female sexual debut by one year for 30% of females that stated they first started having sex at marriage (sensitivity analysis)"))
##         fsd.mar.ind <- which(dat$tfs==dat$tmar) #which sexual debuts occurred at marriage
##         alter.ind <- sample(fsd.mar.ind, size = round(length(fsd.mar.ind)*.3)) # ones to alter
##         dat$tfs[alter.ind] <- dat$tfs[alter.ind] - 12
##     }
##     print("analyzing pooled west african couples data using their respective epidemic curves")
##     save(dat, file = file.path(ndirnm, "wa.Rdata")) # just to indicate which ds working on when looking in file dir
## }                             # if not west africa

## ## Get before couple duration bd, where bd = max(mbd,fbd)
## dat$bd <- apply(cbind(dat$tmar-dat$tms,dat$tmar-dat$tfs), 1, max)
## dat$cd <- dat$tint - dat$tmar ## Get couple duration 
## K <- nrow(dat) ## number of couples
## testpars <- simpars ## for simulation runs
## ## rescale by beta scalars (from args: bmb.sc, bfb.sc, etc...)
## testpars[hazs] <- testpars[hazs] * rep(sapply(paste(hazs, '.sc',sep=''), get))
## ## set up a proposal distribution, Gaussian wit these 'guessed' std devs.
## sd.props <- c(bmb.sd = .001, bfb.sd = .004, 
##               bme.sd = .0015, bfe.sd = .0015,
##               bmp.sd = .003, lrho.sd = .15)
## sd.name <- paste0(parnames, ".sd") ## proposal distr std dev names
## ## load proposal distributions from file if available
## sigma.found <- sum(grepl("sigma", list.files(ndirnm))) > 0 # is there a sigma file?
## if(sigma.found) {
##     print("loading covar matrix for multinormal sampling from file")
##     if("sigmas.Rdata" %in% list.files(ndirnm))        {
##         load(file.path(ndirnm,"sigmas.Rdata"))              # load from array from old simulations
##         sigma <- sigmas[,,group.ind]        #  choose covar matrix for this data set (from earlier adaptations)
##     }else{
##         if(sigma.ad) {                    # use sigma from previous adaptive phase
##             print("loading covar matrix for multinormal sampling from file from previous adaptive phase")
##             load(file.path(ndirnm,"sigma.ad.Rdata")) # load single sigma from earlier adaptive phase in this folder, it's named sigma.ad
##             sigma <- sigma.ad
##         }else{
##             print("loading covar matrix for multinormal sampling from file from previous adaptive phase")
##             load(file.path(ndirnm,"sigma.sm.Rdata")) # load single sigma from earlier sampling phase in this folder, it's named sigma
##             sigma <- sigma.sm
##         }      
##     }
## }else{
##     if(adapt)   print("no covar matrix found in file folder, starting with block univariate sampling")
##     sigma <- NA
## }

## save.image('prepProf.Rdata')

setwd('/home1/02413/sbellan/DHSProject/DHSFitting/')
load(file = 'prepProf.Rdata')
source('DHSFitFunctions.R') # fitting functions
require(Rcpp)
sourceCpp('gonz.cpp')
## test <- dat[sample(1:nrow(dat),10^4, repl=T),]
## test <- dat


sourceCpp('gonz.cpp')
source('DHSFitFunctions.R') # fitting functions

ncpls <- 3000
test <- dat[sample(1:nrow(dat),ncpls, repl=T),]
pre.prepout <- pre.prep(test)
within.prepout <- within.prep(test)
for(nm in names(pre.prepout)) assign(nm, pre.prepout[[nm]]) ## make each of these global for easier access
for(nm in names(within.prepout)) assign(nm, within.prepout[[nm]])

tr <- system.time( cur <- pcalc(simpars, test, browse=F,uncond.mort = T, keep.seros=T, acute.sc = 7))
tc <- system.time(cuc <- pcalcC(max_bd =  max(test$bd), max_cd = max(test$cd), datser = test$ser,
                           pre_fprev = pre.fprev, pre_mprev = pre.mprev, within_fprev = within.fprev, within_mprev = within.mprev,
                           pre_msurv = pre.msurv, pre_fsurv = pre.fsurv, within_msurv = within.msurv, within_fsurv = within.fsurv,
                           within_art_cov = within.art.cov,
                           PreActiveList = PreActiveList, WithinActiveList = WithinActiveList,
                           bmb = simpars['bmb'], bfb = simpars['bfb'],
                           bme = simpars['bme'], bfe = simpars['bfe'],
                           bmp = simpars['bmp'], lrho = simpars['lrho'], lrho_sd = 1/2, trans_ratio = 1,
                           acute_sc = 7, partner_arv = partner.arv, cov_scalar = 1))
colnames(cuc$seros) <- colnames(cur$seros)
identical(cur$seros,cuc$seros)
rbind(tr, tc, tr/tc)
c(R = cur$lprob,C = cuc$lprob)

nits <- 10
new.sampR <- sampler(sd.props=sd.props, dat = test, inits = simpars[-7], niter = nits, nburn = 1, br=F, keep.seros = T, acute.sc = 7, uncond.mort = T, useRcpp=F)
new.sampC <- sampler(sd.props=sd.props, dat = test, inits = simpars[-7], niter = nits, nburn = 1, br=F, keep.seros = T, acute.sc = 7, uncond.mort = T, useRcpp=T)
identical(new.sampR,new.sampC)


sourceCpp('gonz.cpp')
source('DHSFitFunctions.R') # fitting functions

copyv(3)

blah <- matrix(rnorm(100),10,10)
blah
copyc(blah)
identical(blah, copyc(blah))

numcouples <- 10^c(2:5)
times <- numeric(length(numcouples))
for(cc in 1:length(numcouples)) {
    nn <- numcouples[cc]
    maxT <- 800
    examplist <- list(NA)
    for(ii in 1:maxT) examplist[[ii]] <- sample(0:(nn-1), max(101-ii,1))
    pre_fprev <- matrix(runif(nn*maxT),nn,maxT)
    pre_msurv <- matrix(runif(nn*maxT),nn,maxT)
    bmb <-  .05
    times[cc] <- system.time(cplC(max_bd = maxT, pre_fprev = pre_fprev, pre_msurv = pre_msurv, PreActiveList = examplist, bmb = bmb))[3]
}
cbind(numcouples,times, times/numcouples)
