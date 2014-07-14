####################################################################################################
## Analysis of DHS Couple HIV Data
####################################################################################################
## Steve Bellan, September 2012-2013
## steve.bellan@gmail.com
####################################################################################################
## The goal of this analysis is to determine probability of males & females being infected from
## within and outside their relationships as well as estimate transmission coefficients for pre-,
## extra-, & within-couple HIV transmission for both genders for each country.
####################################################################################################
## This script is called by R CMD BATCH with a set of arguments. See DHSFitMK.R which makes
## DHSFitControlFile.txt, the latter has one line with a set of arguments to be sent to a cluster.
####################################################################################################
rm(list=ls())

## setwd('/home1/02413/sbellan/DHSProject/DHSFitting/')
args <- commandArgs(TRUE) ## extract R CMD BATCH arguments
if(length(args)>0) { for(ii in 1:length(args))  eval(parse(text=args[[ii]])) }
source('../SDPSimulations/SimulationFunctions.R') # simulating functions (from other project folder)
source('AuxFxns.R') # other plotting, output manipulating functions
source('DHSFitFunctions.R') # fitting functions
## Load various data inputs (couple copulas, DHS/AIS data, country names, previously fitted parameters, HIV survival model, epidemic curves)
data.fls <- paste0('data files/', c('copula sigmas','allDHSAIS','ds.nm.all','pars.arr.ac','csurv','epic'), '.Rdata')
lapply(data.fls,load,.GlobalEnv)
hazs <- c("bmb","bfb","bme","bfe","bmp","bfp") ## transmission coefficient names
parnames <- hazs; parnames[6] <- 'lrho' # parameters to fit (fitting M->F/F->M transmission ratio, rho)
trans.ratio <- 1 ## m->f : f<-m; within couples. This is the geometric mean of the prior
lrho.sd <- 1/2
print(batchdirnm) ## show where are we saving things.

## These parameters will be used for simulating data that are then fit.
simpars <- c(bmb = .01, bfb = .03, bme = .012, bfe =.01, bmp=.03, lrho = 0)
simpars <- c(simpars, bfp = as.numeric(simpars['bmp'] * exp(simpars['lrho']))) # bfp as a function of bmp & rho
if(simul) { ## if simulating
    ## Set up directory for this simulation
    group <- levels(dat$group)[group.ind]
    dirnm <- paste0(group,'-Sim')
    step <- 1
    ndirnm <- paste0(dirnm,'-', step)
    while(file.exists(ndirnm)) {    ## don't overwrite old files
        ndirnm <- paste0(dirnm,'-', step)
        step <- step + 1
    }
    ndirnm <- file.path(batchdirnm, ndirnm)
    if(!file.exists(ndirnm)) dir.create(ndirnm) ## make directory
    ## Run simulation (see 'SDPSimulations/SimulationFunctions.R' for details on psrun())
    sim <- psrun(country = group.ind, infl.fac = infl.fac, maxN = maxN, vfreq = vfreq, 
                 seed = seed.bump,  # seed to set
                 s.demog = s.demog,        # what country to use for demograhy
                 pars = simpars,
                 psNonPar = psNonPar,
                 last.int = F, sample.tmar = sample.tmar, # non-parametric approach
                 tmar = tmar, each = each, tint = tint, # parametric ps-pop (copulas)
                 death = death,                         # include death?
                 acute.sc = acute.sc,    # acute phase
                 late.sc = late.sc,
                 aids.sc = aids.sc,
                 bmb.sc = bmb.sc, bfb.sc = bfb.sc,      # scale of premarital transmission parameters
                 bme.sc = bme.sc, bfe.sc = bfe.sc,      # scale of extramarital transmission parameters
                 bmp.sc = bmp.sc, bfp.sc = bfp.sc, #  within-couple transmission coefficients, to scale fitted ones                  
                 het.b = het.b, het.b.sd = het.b.sd, het.b.cor = het.b.cor,
                 het.e = het.e, het.e.sd = het.e.sd, het.e.cor = het.e.cor,
                 het.p = het.p, het.p.sd = het.p.sd, het.p.cor = het.p.cor, 
                 het.gen = het.gen, # genetic heterogeneity
                 het.gen.sd = het.gen.sd,
                 het.gen.cor = het.gen.cor,
                 het.beh = het.beh, # behavioral heterogeneity (only applies to premarital & extramarital)
                 het.beh.sd = het.beh.sd,
                 het.beh.cor = het.beh.cor,
                 scale.by.sd = scale.by.sd, # adjust beta means to keep geometric mean constant with increasing het
                 scale.adj = scale.adj,   # adjust them arbitrarily
                 out.dir = ndirnm,  nc = nc, make.jpgs = F,
                 browse = F) 
    load(sim) ## sim is a file name, load it
    dat <- output$evout ## store line list as dat (to be fit)
    dat <- add.cat(dat) ## add some more variables to it
    simdat <- dat
    dat <- dat[dat$alive,]  ## only fit models to those who lived to end of survey
}else{ ## otherwise fitting real data
    ## Set up new directory to store results
    group <- levels(dat$group)[group.ind]
    dirnm <- paste0(group,'RealAc',acute.sc, '-jb', jobnum)
    step <- 1
    ndirnm <- paste0(dirnm,'-', step)
    while(new.folder & file.exists(ndirnm)) { ## don't continue from old analysis (otherwise redo new folder)
        ndirnm <- paste0(dirnm,'-', step)
        step <- step + 1
    }
    ndirnm <- file.path(batchdirnm, ndirnm)
    if(!file.exists(ndirnm))        dir.create(ndirnm) ## create directory
    wa <- levels(dat$group)[group.ind]=="WA"
    if(!wa) { ## if not WA just set countries
        group <- levels(dat$group)[group.ind]
        print(paste(nlevels(dat$ds),"country data sets"))
        ds <- levels(dat$ds)[grepl(group, levels(dat$ds))]   # data set we're working on
        dat <- dat[dat$ds %in% ds,]           # only looking at that data set
        print(paste("analyzing",paste(ds, collapse = " & "),"couples data using", dat$epic.nm[1],
                    "epidemic curve"))
        if(fsd.sens) { # if doing sensitivity analysis to female sexual debut
            ## Of females saying their sexual debut occured at
            ## marriage, assume 30% were lying and that actually
            ## sexual debut occurred one year earlier.
            print(paste("lowering female sexual debut by one year for 30% of females that stated they first started having sex at marriage (sensitivity analysis)"))
            fsd.mar.ind <- which(dat$tfs==dat$tmar) #which sexual debuts occurred at marriage
            alter.ind <- sample(fsd.mar.ind, size = round(length(fsd.mar.ind)*.3)) # ones to alter
            dat$tfs[alter.ind] <- dat$tfs[alter.ind] - 12
        }
        save(dat, file = file.path(ndirnm, paste0(group,".Rdata"))) # just to indicate which ds working on when looking in file dir
    }else{ ## For West Africa, need to select all countries in region
        group <- levels(dat$group)[group.ind]
        print(paste(nlevels(dat$ds),"country data sets"))
        dat <- dat[dat$group == group,]                 # west africa
        if(fsd.sens) { # if doing sensitivity analysis to female sexual debut
            ## Of females saying their sexual debut occured at
            ## marriage, assume 30% were lying and that actually
            ## sexual debut occurred one year earlier.
            print(paste("lowering female sexual debut by one year for 30% of females that stated they first started having sex at marriage (sensitivity analysis)"))
            fsd.mar.ind <- which(dat$tfs==dat$tmar) #which sexual debuts occurred at marriage
            alter.ind <- sample(fsd.mar.ind, size = round(length(fsd.mar.ind)*.3)) # ones to alter
            dat$tfs[alter.ind] <- dat$tfs[alter.ind] - 12
        }
        print("analyzing pooled west african couples data using their respective epidemic curves")
        save(dat, file = file.path(ndirnm, "wa.Rdata")) # just to indicate which ds working on when looking in file dir
    }                             # if not west africa
}                                 # end if/else simulation statement
save.image(file=file.path(ndirnm,"pars workspace.Rdata"))

## Get before couple duration bd, where bd = max(mbd,fbd)
dat$bd <- apply(cbind(dat$tmar-dat$tms,dat$tmar-dat$tfs), 1, max)
dat$cd <- dat$tint - dat$tmar ## Get couple duration 
K <- nrow(dat) ## number of couples
testpars <- simpars ## for simulation runs
## rescale by beta scalars (from args: bmb.sc, bfb.sc, etc...)
testpars[hazs] <- testpars[hazs] * rep(sapply(paste(hazs, '.sc',sep=''), get))
## set up a proposal distribution, Gaussian wit these 'guessed' std devs.
sd.props <- c(bmb.sd = .001, bfb.sd = .004, 
              bme.sd = .0015, bfe.sd = .0015,
              bmp.sd = .003, lrho.sd = .15)
sd.name <- paste0(parnames, ".sd") ## proposal distr std dev names
## load proposal distributions from file if available
sigma.found <- sum(grepl("sigma", list.files(ndirnm))) > 0 # is there a sigma file?
if(sigma.found) {
    print("loading covar matrix for multinormal sampling from file")
    if("sigmas.Rdata" %in% list.files(ndirnm))        {
        load(file.path(ndirnm,"sigmas.Rdata"))              # load from array from old simulations
        sigma <- sigmas[,,group.ind]        #  choose covar matrix for this data set (from earlier adaptations)
    }else{
        if(sigma.ad) {                    # use sigma from previous adaptive phase
            print("loading covar matrix for multinormal sampling from file from previous adaptive phase")
            load(file.path(ndirnm,"sigma.ad.Rdata")) # load single sigma from earlier adaptive phase in this folder, it's named sigma.ad
            sigma <- sigma.ad
        }else{
            print("loading covar matrix for multinormal sampling from file from previous adaptive phase")
            load(file.path(ndirnm,"sigma.sm.Rdata")) # load single sigma from earlier sampling phase in this folder, it's named sigma
            sigma <- sigma.sm
        }      
    }
}else{
    if(adapt)   print("no covar matrix found in file folder, starting with block univariate sampling")
    sigma <- NA
}

## Prepare data set for fitting. 
pre.prepout <- pre.prep(dat)
within.prepout <- within.prep(dat)
for(nm in names(pre.prepout)) assign(nm, pre.prepout[[nm]]) 
for(nm in names(within.prepout)) assign(nm, within.prepout[[nm]])

start.time1 <- Sys.time() ## start time of sampling
if(all.cores) {           ## make wrapper function around sampler for mclapply
    wrp <- function(seed=1, multiv=F, covar=NULL, acute.sc, dat, lrho.sd,
                    niter, survive, browse,
                    nthin, nburn)
      { ## new version of wrp needs to save progress for long chains (WA)
        inits.temp <- init.fxn(seed = seed) # initial conditions different for each seed
        sampler(dat = dat, sd.props = sd.props, inits = inits.temp, acute.sc = acute.sc, browse = browse,
                multiv = multiv, covar = covar, lrho.sd = lrho.sd,
                verbose = T, tell = tell, seed = seed,
                niter = niter, survive = survive,
                nthin = nthin, keep.seros = TRUE,
                nburn = nburn)
      }
    if(adapt) { ## Adaptive (d.) phase to get multivariate normal sampler
        print("beginning adaptive phase")
        d.out <- mclapply((seed.bump + 1:nc), wrp, 
                          dat = dat, acute.sc = acute.sc, multiv = sigma.found, covar = sigma, browse=F, lrho.sd = lrho.sd,
                          survive = survive, niter = d.niter, nthin = d.nthin, nburn = d.nburn)
        save.image(file=file.path(ndirnm,"workspace.Rdata"))
        ## reformat into mcmc object
        mcmc.d.out <- list(NA)
        d.aratio <- 0
        for(ii in 1:nc) {
            mcmc.d.out[[ii]] <- as.mcmc(d.out[[ii]]$out)
            d.aratio <- d.aratio + d.out[[ii]]$aratio
            if(ii==1) { ## initialize
                init.adapt <- d.out[[ii]]$inits
            }else{ ## append
                init.adapt <- rbind(init.adapt, d.out[[ii]]$inits)
            }
        }
        mcmc.d.out <- mcmc.list(mcmc.d.out) # as mcmc.list
        d.aratio <- d.aratio/nc             # acceptance ratio
        print(paste("adaptive phase aratio is", round(d.aratio,2)))
        for(pp in parnames) assign(paste0(pp,'.vec'), unlist(mcmc.d.out[,pp])) ## pull out all chains into vectors
        bfp.vec <- bmp.vec * exp(lrho.vec)  
        posts <- data.frame(bmb = bmb.vec, bfb = bfb.vec, bme = bme.vec, bfe = bfe.vec, bmp = bmp.vec, lrho = lrho.vec)
        sbpairs(posts, truepars = testpars[parnames], show.lines = simul, ## plot posterior correlations after adaptive phase
                file.nm = file.path(ndirnm,"posterior pairs after adaptive phase"), width = 12, height = 12,
                cex = 1, col = "black", nrpoints = 200, do.jpeg = T)
        mu <- colMeans(posts) ## posterior means
        sigma.ad <- cov.wt(posts)$cov ## posterior covariance matrix, then plot what proposal distr is gonna look like
        sbpairs(rmnorm(3000, mean = mu, varcov = sigma.ad), truepars = testpars[parnames], show.lines = simul,
                file.nm = file.path(ndirnm, "adapted proposal distr after adaptive phase"), width = 12, height = 12,
                cex = 1, col = "black", nrpoints = 200, do.jpeg = T, do.pdf = F)
        save.image(file=file.path(ndirnm,"workspace.Rdata")) ## backup workspace
        save(sigma.ad, file=file.path(ndirnm,"sigma.ad.Rdata")) ## save covariance matrix from adaptive phase
        sigma <- sigma.ad ## update sigma to adapted sigma
      } # end adaptive phase, sigma has been updated if this is run, otherwise it's been loaded
    print("beginning sampling")
    ## use multicore to run a chain on each core with seeds 1:nc (where nc is number chains/core)
    out <- mclapply((seed.bump + 1:nc), wrp, acute.sc = acute.sc, multiv = T, covar = sigma, browse = F,
                    survive = survive, niter = niter, nthin = nthin, nburn = nburn) ## this takes a long time!
    save.image(file=file.path(ndirnm,"workspace.Rdata"))  ## backup workspace
    mcmc.out <- list(NA)   ## reformat into mcmc object
    aratio <- 0
    for(ii in 1:nc) { ## combine each core's chain into mcmc.list
        mcmc.out[[ii]] <- as.mcmc(t(out[[ii]][[1]]))
        aratio <- aratio + out[[ii]]$aratio
        if(ii==1) {
            init.samp <- out[[ii]]$inits
        }else{
            init.samp <- rbind(init.samp, out[[ii]]$inits)
        }
    }
    mcmc.out <- mcmc.list(mcmc.out)
    aratio <- aratio/nc
    for(pp in parnames) assign(paste0(pp,'.vec'), unlist(mcmc.out[,pp])) ## pull out all chains into vectors
    bfp.vec <- bmp.vec * exp(lrho.vec)  
    posts <- data.frame(bmb = bmb.vec, bfb = bfb.vec, bme = bme.vec, bfe = bfe.vec, bmp = bmp.vec, lrho = lrho.vec)
    sigma.sm <- cov.wt(posts)$cov #estimate covariance matrix, then plot what proposal distr is gonna look like
    save(sigma.sm, file = file.path(ndirnm,"sigma.sm.Rdata"))
    sbpairs(posts, truepars = testpars[parnames], show.lines = simul, file.nm = file.path(ndirnm,"pairs after sample phase"),
            width = 12, height = 12, cex = 1, col = "black", nrpoints = 200, do.jpeg = T)
  }else{                                # if just doing one core
    inits <- init.fxn(seed = seed)
    out <- sampler(sd.props = sd.props, inits = inits, acute.sc = acute.sc,
            verbose = T, tell = 20, seed = seed, lrho.sd = lrho.sd,
            niter = niter, survive = survive,
            nthin = nthin,
            nburn = nburn)
    mcmc.out <- as.mcmc(t(out[[1]]))
    aratio <- out$aratio
  }  
save.image(file=file.path(ndirnm,"workspace.Rdata"))
save(aratio, file = file.path(ndirnm,paste("aratio is", round(aratio,2), ".Rdata")))
print(paste("aratio is", signif(aratio,2)))
####################################################################################################
## END MCMC SAMPLING

###################################################################### 
## Extract parameters' posterior and save to file (use parallel processing)
sum.wrap <- function(col.ind, xx) {
    summary(xx[,col.ind], quantiles = c(.025, .5, .975))$quant
  }
pars <- abind(mclapply(1:ncol(mcmc.out[[1]]), sum.wrap,  xx = mcmc.out), along = 0, new.names=colnames(mcmc.out[[1]]))
## never save over old files
file.name <- paste("pars", format(Sys.time(), "%Y%m%d"), sep = "-")
ii <- 1
while(file.exists(file.path(ndirnm,paste0(file.name, ".Rdata")))) {
    file.name <- paste(file.name, "-",ii)
    ii <- ii+1
}
save(pars, file = file.path(ndirnm,paste0(file.name, ".Rdata")))
save(mcmc.out, file = file.path(ndirnm,paste0(file.name, "chains.Rdata")))
out.csv <- signif(pars,3)
out.csv <- data.frame(median = out.csv[,2], CI95 = paste("(",out.csv[,1],", ",out.csv[,3],")", sep=""))
write.csv(out.csv, file=file.path(ndirnm,"out.csv"))
save.image(file=file.path(ndirnm,"workspace.Rdata"))

####################################################################################################
####################################################################################################
## ALL CODE BELOW IS FOR FIGURES
####################################################################################################
####################################################################################################
###################################################################### 
## Figure 0 - MCMC diagnostics (for beta's only)
pdf(file.path(ndirnm, "Fig 0 - mcmc diagnostics.pdf"))
show.cols <- colnames(mcmc.out[[1]]) %in% c(parnames,"bfp")
plot(mcmc.out[,show.cols])
dev.off()
######################################################################
## Gelman-Rubin diagnotics (save to file)
library(coda)
betpars <- c("bmb", "bfb", "bme", "bfe", "bmp", "bfp")
beta.out <- mcmc.out[, betpars, drop=FALSE]
gelout <- gelman.diag(beta.out)
gelout
###################################################################### 
## Figure 1 - plot probability of transmission by route over time
pdf(file.path(ndirnm, "Fig 1 - risk by route.pdf"), width = 7, height = 8)
par(mfrow=c(3,2))
nn <- nrow(dat)
xs <- 1:40 ## years of sex
xlim <- c(0, 46)
xlabs <- rep(c("years sexually active before relationship", rep("relationship duration",2)), 2)
mains <- c(paste(c("male","female"), "risk of infection prior to relationship"),
          paste(c("male","female"), "risk of infection from extra-couple sex"),
          paste(c("male","female"), "risk of infection from within-couple sex"))
for(hh.i in 1:6) {
    hh <- hazs[hh.i]
    plot(0,0, xlim = xlim, ylim = c(0,1),
         type = "n", bty = "n",
         ylab = expression(1-exp(-beta[M]*(x[M,i]-r[i]))), xlab = xlabs[hh.i],
         main = mains[hh.i])
    polygon(c(xs, rev(xs)),
            c(1-exp(-xs*pars[hh,"2.5%"]), rev(1-exp(-xs*pars[hh,"97.5%"]))),
            col = "gray", border = NA)
    lines(xs, 1-exp(-xs*pars[hh,"50%"]), lwd = 2)
}
dev.off()

###################################################################### 
## Figure 2 - prob an infection is from extracouple sex by serostatus
## by ysa and reldur
###################################################################### 
cex <- .8
medpars <- pars[parnames,2]
pis <- pcalc(medpars, dat = dat, trace = T, give.pis=T, survive=survive, lrho.sd = log(trans.ratio))$pis
######################################################################
breaks <- seq(0,1, by = .1)
xlim <- c(0, 35)
ylim <- c(0, 35)
xlab <- "years sexually active before relationship"
ylab <- "relationship duration"
mains <- c("M+F+","M+F-","M-F+","M-F-")
######################################################################
cex <- .65
rmp <- colorRamp(c("yellow","red"))     #create color ramp
pal <- colorRampPalette(c("yellow","red"))     #create color palette (for legend)
## pdf(file.path(ndirnm, "Fig 2 - serostatus by ysa prior and reldur normalized for prev with epidemic curve.pdf"),
##     width = 6.5, height = 4
tiff <- F
if(tiff) { ## use tiff to keep file size small for lots of points (& for publication)
    tiff(file.path(ndirnm,"Fig 2 - serostatus by ysa prior and reldur normalized for prev with epidemic curve.tiff"),
         width = 5.5, height = 4, units = "in", res = 300)
  }else{
    pdf(file.path(ndirnm, "Fig 2 - serostatus by ysa prior and reldur normalized for prev with epidemic curve.pdf"),
         width = 5.5, height = 4)
  }
## do it for each data set since the interview times were different for same countries and so the
## plot shold show how long the couples were together too
par(mar = c(.5,1,1,.5), oma = c(5,5,0,9))
for(dd in unique(dat$epic.nm)) {
    cc <- dat$epic.ind[dat$epic.nm==dd][1]   #find epidemic curve for that data set
    epic.col <- "blue"
    layout(t(matrix(1:4,2,2)))
    ylim <- c(0, 30) #max(dat$m.bef.pm, dat$f.bef.pm))/12
    xlim <- c(1975,2012)
    ylab <-  "YSA before couple formation"
    xlab <- "date of couple formation"
    mains <- c("M+F+","M+F-","M-F+","M-F-")
    mains <- c("B","A","","C")
    for(ii in 2:1) { ## for mSDC & ++ serostatus couples
        if(ii!=4) cols.show <- rgb(rmp(pis$piCe.A[dat$ser==ii & dat$epic.ind==cc]), max = 255)
        if(ii==4) cols.show <- "black"
        plot(1900 + 1/12*(dat$tint-dat$mardur.mon)[dat$ser==ii & dat$epic.ind==cc],
             1/12*(dat$tmar-dat$tms)[dat$ser==ii & dat$epic.ind==cc],
             col = cols.show, xlab = "",
             pch = 19, cex = cex, axes = F,
             xlim = xlim, ylim = ylim, bty = "n",
             main = mains[ii])
        xs <- epicm[,1]/12 + 1900
        if(ii==2) axis(2, at = seq(0, 30, by = 10), las = 2)
        if(ii==1) axis(2, at = seq(0, 30, by = 10), labels = NA, las = 2)
        axis(1, at = seq(1980, 2010, by = 10), labels = NA)
        lines(xs[xs>1975], epicf[xs>1975, cc]*max(ylim), col = epic.col, lwd = 1)
        if(ii==1) axis(4, at = seq(0, max(ylim), l=5), seq(0, 100, l = 5), col = epic.col, las = 2)
      }
    mains <- c("F+M+","","F+M-","F-M-")
    mains <- c("D","","C","F")    
    for(ii in c(3,1)) { ## for fSDC & ++ serostatus couples
        if(ii!=4) cols.show <- rgb(rmp(pis$piC.eA[dat$ser==ii & dat$epic.ind==cc]), max = 255)
        if(ii==4) cols.show <- "black"
        plot(1900 + 1/12*(dat$tint-dat$mardur.mon)[dat$ser==ii & dat$epic.ind==cc],
             1/12*(dat$tmar-dat$tfs)[dat$ser==ii & dat$epic.ind==cc],
             ##          1/12*(dat$f.bef.pm)[dat$ser==ii],
             col = cols.show, pch = 19, las = 2, axes = F,
             xlim = xlim, ylim = ylim, bty = "n", xlab = "", cex = cex,
             main = mains[ii])
        if(ii==3) axis(2, at = seq(0, 30, by = 10), las = 2)
        if(ii==1) axis(2, at = seq(0, 30, by = 10), labels = NA, las = 2)
        axis(1, at = seq(1980, 2010, by = 10), las = 2)
        lines(xs[xs>1975], epicm[xs>1975, cc]*max(ylim), col = epic.col, lwd = 1)
        if(ii==1) axis(4, at = seq(0, max(ylim), l=5), seq(0, 100, l = 5), col = epic.col, las = 2)
      }
    cex.ax <- .7
##     mtext("male YSA before couple formation", side = 2, outer = T, line = 2, cex = cex.ax, adj = .92)
##     mtext("female YSA before couple formation", side = 2, outer = T, line = 2, cex = cex.ax, adj = .08)
##     mtext(paste(xlab), side = 1, outer = T, line = 2, cex = cex.ax)
##     mtext("F HIV pop. prevalence", side = 4, outer = T , line = 2, adj = .9, cex = cex.ax)
##     mtext("M HIV pop. prevalence", side = 4, outer = T , line = 2, adj = .15, cex = cex.ax)
##     mtext(dd, side = 3, outer = T, line = .5) # data set title
    ##       mtext(paste("male",xlab), side = 1, outer = T, line = -21, cex = 1.2)
  } # end loop through pooled data sets
dev.off()
## make tiff legend separate for combination in Illustrator
tiff(file.path(ndirnm,"fig 2 legend.tiff"), width = .9, height = 3, units = "in", res = 300)
par(mar=rep(0,4))
plot(0,0,type="n",axes=F, xlim = c(-.1,.2), ylim = c(-.1,.9))
cols <- pal(100)
color.legend(0,.1,.05,.8, seq(0,1, length=11), rect.col = cols, gradient = "y", cex = .8)
dev.off()
## PDF legend
pdf(file.path(ndirnm, "fig 2 legend.pdf"), width = .9, height = 3)
par(mar=rep(0,4))
plot(0,0,type="n",axes=F, xlim = c(-.1,.2), ylim = c(-.1,.9))
cols <- pal(100)
color.legend(0,.1,.05,.8, seq(0,1, length=11), rect.col = cols, gradient = "y", cex = .8)
dev.off()
######################################################################

###################################################################### 
## Figure 3 - Hazard posteriors, medians & 95% CIs
###################################################################### 
pdf(file.path(ndirnm, "Fig 3 - hazard posteriors.pdf"), width = 6, height = 3.5)
par(mar=c(5,10,0.5,.5))
show <- match(hazs, rownames(pars))
labs <- c("M before relationship", "F before relationship",
          "M extra-couple sex", "F extra-couple sex",
          "M from partner", "F from partner")
plot(12*pars[show,2], 6:1, pch = 15, cex = 1,
     ylim = c(.8,6.2), xlim=c(0,12*max(pars[show,])), bty = "n", yaxt = "n",
     ylab="", xlab = expression(beta[yearly]))
arrows(12*pars[show,1],6:1,  12*pars[show,3],  6:1, angle = 90, length=.1, code = 3, lwd = 2)
axis(2, at = 6:1, label = labs, las = 2, cex=2)
## axis(1, at = seq(0,.3, by = .05))
dev.off()
###################################################################### 

###################################################################### 
## Figure 6 - log(bfp/bmp) posterior & prior
###################################################################### 
pdf(file.path(ndirnm, "Figure 6 - lmf to lfm post.pdf"), width = 6, height = 5)
par(mar = c(4.5,4,2.5,0))
xx <- seq(0, 10, length.out=1000)
yy <- dlnorm(xx, mean = log(trans.ratio), sd = 1/2)
hist(exp(unlist(mcmc.out[,"lrho"])), breaks = 100, col="gray", border = NA,
     xlab = expression(beta[mf]/beta[fm]), ylab = "probability density", main = "",
     xlim = c(0,10), freq = F, ylim = c(0, 1))
lines(xx, yy, lwd = 2)
legend("topright", leg = c("prior","posterior"), col = c("black","gray"), lwd = 2, bty = "n")
dev.off()
###################################################################### 
## Figure 7 - bfp/bmp prior (only
###################################################################### 
pdf(file.path(ndirnm, "Figure 7 - prior beta out.pdf"), width =4, height = 3.5)
par(mar = c(4.5,4,2.5,0))
xx <- seq(0, 6, length.out=1000)
yy <- dlnorm(xx, mean = log(trans.ratio), sd = 1/2)
plot(xx,yy, lwd = 2, type = "l",
     xlab = expression(beta[Fpartner]/beta[Mpartner]), ylab = "probability density", main = "",
     xlim = c(0,6), ylim = c(0, 1), bty = "n")
dev.off()
###################################################################### 
## Figure 8 - years sexally active before relationship distribution
###################################################################### 
pdf(file.path(ndirnm, "Figure 8 - ysa for disc couples.pdf"), width = 6, height = 8)
par(mfrow=c(2,1))
hist(1/12*(dat$tmar-dat$tms)[dat$ser==2], breaks = seq(0,55,by=1/2), xlab = "MYSA before relationship",
     main="M+ F- couples", col = "black")
hist(1/12*(dat$tmar-dat$tfs)[dat$ser==3], breaks = seq(0,55,by=1/2), xlab = "FYSA before relationship",
     main="M- F+ couples", col = "black")
dev.off()
###################################################################### 

###################################################################### 
## Figure 9 - beta ratios
###################################################################### 
pdf(file.path(ndirnm, "Figure 9 - beta ratios.pdf"), width = 5, height = 3.5)
par(mar=c(4,7,.5,.5))
show <- match(c("rr.mf.bef", "rr.mf.exc",
                "rr.m.eb", "rr.f.eb"), rownames(pars)) # SHOWING THE INVERSES (SEE LABELS)
labs <- c(expression(beta[Fbefore] / beta[Mbefore]),
          expression(beta[Fextra] / beta[Mextra]),
          expression(beta[Mbefore] / beta[Mextra]),
          expression(beta[Fbefore] / beta[Fextra]))
plot(1/pars[show,2], 4:1, pch = 15, cex = 2, log="x",
     ylim = c(.8,4.2), xlim=c(.05,20), bty = "n", axes=F, ylab="", xlab = "rate ratio")
segments(1, 4.2, 1, .8, lty = 2, lwd = 2)
arrows(1/pars[show,1],4:1,  1/pars[show,3],  4:1, angle = 90, length=.1, code = 3, lwd = 2)
axis(2, at = 4:1, label = labs, las = 2, cex=2)
axis(1, at = c(.05,.1,.2,.5,1,2,5,10,20), las = 2,
     label = c("1/20","1/10","1/5","1/2","1","2","5","10","20"))
dev.off()
######################################################################
###################################################################### 
## Figure 9b - beta ratios
###################################################################### 
pdf(file.path(ndirnm, "Figure 9b - contact ratios.pdf"), width = 5, height = 3.5)
par(mar=c(4,7,.5,.5))
show <- match(c("rr.mf.bef.cont", "rr.mf.exc.cont",
                "rr.m.out", "rr.f.out"), rownames(pars)) # SHOWING THE INVERSES (SEE LABELS)
labs <- c(expression(beta[Fbefore] / beta[Mbefore]),
          expression(beta[Fduring] / beta[Mduring]),
          expression(beta[Mbefore] / beta[Mduring]),
          expression(beta[Fbefore] / beta[Fduring]))
plot(1/pars[show,2], 4:1, pch = 15, cex = 2, log="x",
     ylim = c(.8,4.2), xlim=c(.05,20), bty = "n", axes=F, ylab="", xlab = "rate ratio")
segments(1, 4.2, 1, .8, lty = 2, lwd = 2)
arrows(1/pars[show,1],4:1,  1/pars[show,3],  4:1, angle = 90, length=.1, code = 3, lwd = 2)
axis(2, at = 4:1, label = labs, las = 2, cex=2)
axis(1, at = c(.05,.1,.2,.5,1,2,5,10,20), las = 2,
     label = c("1/20","1/10","1/5","1/2","1","2","5","10","20"))
dev.off()
######################################################################

###################################################################### 
pdf(file.path(ndirnm, "Figure 11 - # next year of transmission.pdf"), width = 4, height = 5)
  ## Calculate probability each uninfected person is infected in next
  ## year using estimated beta's as well as current prevalence.
  par(mar=c(10,4,0.5,.5))
show <- match(c("n.m.part.tot", "n.m.exc.tot", "n.f.part.tot", "n.f.exc.tot"),
              rownames(pars))                
nms <- c("male: partner", "male: extracouple", "female: partner", "female: extracouple")
bp <- barplot(pars[show,2], names.arg = nms, col = c("blue","red","blue","red"),
ylab = "# incident infections in next 12 months",
las = 2, ylim = c(0, max(pars[show,])))
arrows(bp,pars[show,1],bp,pars[show,3], angle = 90, length=.1, code = 3, lwd = 2)
mtext(paste("N =",c(sum(dat$ser %in% c(3:4)))), side = 1, adj = .22, line = 0)
mtext(paste("N =",c(sum(dat$ser %in% c(2,4)))), side = 1, adj = .82, line = 0)
dev.off()
######################################################################

###################################################################### 
pdf(file.path(ndirnm, "Figure acute - proportion of within transm due to acute.pdf"), width = 4, height = 5)
  ## Calculate probability each uninfected person is infected in next
  ## year using estimated beta's as well as current prevalence.
  par(mar=c(10,4,0.5,.5))
show <- match(c("pipUaA","pipUcA", "piUpaA","piUpcA"), rownames(pars))                
nms <- c("M acute from partner", "M chronic from partner", "F acute from partner", "F chronic from partner")
bp <- barplot(pars[show,2], names.arg = nms, col = c("blue","dark blue","blue","dark blue"),
ylab = "proportion of all observed infections",
las = 2, ylim = c(0, max(pars[show,])))
arrows(bp,pars[show,1],bp,pars[show,3], angle = 90, length=.1, code = 3, lwd = 2)
dev.off()
######################################################################
    
###################################################################### 
pdf(file.path(ndirnm, "Figure 12 - prob any new inf is extracouple.pdf"), width = 3, height = 5)
## Calculate probability each uninfected person is infected in next
## year using estimated beta's as well as current prevalence.
par(mar=c(3,4,0.5,.5))
nms <- c("male","female")
show <- match(c("prop.exc.m", "prop.exc.f"), rownames(pars))
bp <- barplot(pars[show,2], names.arg = nms,
        ylab = "probability new infection is from extracouple intercourse",
        las = 1, ylim = c(0, 1))
arrows(bp,pars[show,1],bp,pars[show,3], angle = 90, length=.1, code = 3, lwd = 2)
dev.off()
######################################################################

###################################################################### 
pdf(file.path(ndirnm, "Fig 16 - survival times.pdf"), width = 4.5, height = 3)
par(mar=c(4,4,1,1))
## create discrete cumulative probability of mortality
## it's age dependent, fit by eyeballing to CASCADE study
xseq <- 1:(12*60)
age.seq <- seq(20*12,60*12, by = 10*12)
cols <- c("orange","green","pink","light blue", "dark gray")
for(ii in 1:length(age.seq)) {
    aa <- age.seq[ii]
    shp <- 2.3
   scl <- 2000/shp/(aa/12)^.53
    cmort <- pweibull(xseq, shape = shp, scale = scl)
    csurv.temp <- 1-cmort
    if(ii==1) plot(xseq/12, csurv.temp, type = "l", xlim = c(0,25), col=cols[ii], lwd = 3, yaxt = "n",
         xlab = "years since seroconversion", ylab="probability of survival", bty = "n")
    if(ii!=1) lines(xseq/12, csurv.temp, type = "l", col=cols[ii], lwd = 3)
  }
legend("topright",paste(age.seq/12,"yrs old"), col = cols, pch = 15, bty = "n", cex = 1, title = "age at seroconversion")
axis(2, seq(0,1,l=5), las = 2)
dev.off()
######################################################################

###################################################################### 
ctraj(medpars, dat, browse =F,
      plot.cpls = sample(1:nrow(dat),10,replace=F),
      surv = survive,                   # show surv curves o plot
      nsurv = F,                        # don't show marginal curves 
      lty.surv = 1,                     # line type for surv curves
      dead =T, col.dead = "black", lty.dead = 1,
      survive = survive,                # plot point at survival curve
      pdf.name = file.path(ndirnm, "Figure 14 - Probability trajectories for some couples.pdf"))
###################################################################### 

## compare fitted to actual values from simulated data
if(simul) compsim(simdat, pars, testpars, dirnm = ndirnm, browse = F)

######################################################################
## Print processing time to file
hours <- round(as.numeric(difftime(Sys.time(), start.time1, unit = "hour")),3)
save(hours, file = file.path(ndirnm,paste("took",hours,"hrs.Rdata")))
save.image(file=file.path(ndirnm,"workspace.Rdata")) ## resave workspace
print(getwd())
