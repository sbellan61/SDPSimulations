####################################################################################################
## Makes control files for each analysis within which each line giving one R CMD BATCH command line
## to run on a cluster.
####################################################################################################
#rm(list=ls())                                  # clear workspace
require(data.table)
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/DHSProject/SDPSimulations/')
if(grepl('stevenbellan', Sys.info()['login'])) setwd('~/Documents/R Repos/SDPSimulations/SDPSimulations/')
load("../DHSFitting/data files/ds.nm.all.Rdata") # country names
load('../DHSFitting/data files/pars.arr.ac.Rdata')    # load acute phase relative hazards used to fit (in.arr[,,2])
load('../DHSFitting/data files/CFJobsToDo.Rdata') ## for finishing up jobs from last run that didn't get finished due to cluster problems.
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') #  transmission coefficient names, for convenience
nc <- 12                                       # core per simulation
## source('CounterFactualMK.R')

####################################################################################################
####################################################################################################
##   COUNTERFACTUAL ANALYSIS
####################################################################################################
## Using an HIV transmission model fit to Demographic and Health Surveys in 18 sub-Saharan African
## countries, we used counter-factual simulations to examine how serodiscordant proportions are
## affected by AIDS mortality rates and pre-couple, extra-couple, and within-couple HIV transmission
## rates.
####################################################################################################
countries <- 1:length(ds.nm)

#countries <- which(ds.nm=='Zambia')

each.val <- 200                          #  number of couples per couple formation (marital) cohort
counterf.betas <- F                       # change betas in counterfactuals? if not change beta_within & c's (so beta_within affects all routes)
sub.betas <- F                           # substitute betas? if not beta_within & c's
rtsc <- c(0, 1/10, 1/5, 1/2, 1, 2, 5, 10)  # transmission route scalars
nrtsc <- length(rtsc)                      # how many scalars?
# acutes <- as.numeric(in.arr[,1,2])    # acute phase relative hazards we used to fit in fitting phase
acutes <- c(1,5,7,10)
nac <- length(acutes)                 # how many are there?
hsds <- c(.5,1,2,3)                     # standard deviation of log(hazards)
hsds.rts <- c(0,1,2)                      # same but use smaller subset when also scaling transmission routes
cors <- c(0,.4,.8)                      # inter-partner correlations
ncors <- length(cors)

nn <- 400 # number of simulations per country-acute combination (must be bigger than max(sel) later but is fine to leave big
substitute <- F                         # not a substitution analysis
totn <- 0                               # total number of simulations (steps up to final value)
num.doing <- 0
outdir <- file.path('results','CounterFactual')

nhsds <- length(hsds)
hets <- c('b','e','p','gen','beh')
labs <- c('pre-','extra-','within-', 'genetic', 'behavioral')

currbatch <- 1
blocks <- data.table(batch = 99, simj = 1:160, lab = 'as fitted', death=T, 
                     het.gen=F, het.gen.sd=0, het.gen.cor=0, het.beh=F, het.beh.sd=0, het.beh.cor=0,
                     het.b=F, het.b.sd=0, het.b.cor=0, het.e=F, het.e.sd=0, het.e.cor=0, het.p=F, het.p.sd=0, het.p.cor=0,
                     bmb.sc=1, bfb.sc=1, bme.sc=1, bfe.sc=1, bmp.sc=1, bfp.sc=1,
                     acute.sc=1, late.sc=1, aids.sc=1)
blocks[1, batch:=currbatch]
for(ghsd in hsds.rts) {
    nextsimj <- blocks[batch==currbatch, max(simj)+1]
    currbatch <- currbatch + 1
    blocks[nextsimj, c('batch','lab','death','het.gen', 'het.gen.sd'):=.(currbatch,paste0('no AIDS mortality, gen het=',ghsd), F, ghsd>0, ghsd)]
    for(hi in 1:3) {
        hh <- hets[hi]
        nextsimj <- blocks[batch==currbatch, max(simj)+1]
        currbatch <- currbatch + 1
        blocks[nextsimj:(nextsimj+nrtsc-1), c('batch', 'lab', paste0('bm', hh,'.sc'),paste0('bf', hh,'.sc'), 'het.gen', 'het.gen.sd')
               := .(currbatch, paste0('scale ',labs[hi], '-couple, gen het=', ghsd), rtsc, rtsc, ghsd>0, ghsd)]
    }
}
for(cr in cors) {
    for(hi in 1:5) {
        hh <- hets[hi]
        nextsimj <- blocks[batch==currbatch, max(simj)+1]
        currbatch <- currbatch + 1
        blocks[nextsimj:(nextsimj+nhsds-1), c('batch', 'lab', paste0(paste0('het.',hh), c('','.sd','.cor'))) := .(currbatch, paste0(hh, ' cor=', cr), T, hsds, cr)]
    }
    ## pre-/extra- route heterogeneity (different rfs per route)
    nextsimj <- blocks[batch==currbatch, max(simj)+1]
    currbatch <- currbatch + 1
    blocks[nextsimj:(nextsimj+nhsds-1), c('batch', 'lab', paste0(rep(c('het.b','het.e'), each = 3), c('','.sd','.cor')) ) := .(currbatch, paste0('pre/extra cor=', cr), T, hsds, cr, T, hsds, cr)]
    ## all route heterogeneity (different rfs per route)
    nextsimj <- blocks[batch==currbatch, max(simj)+1]
    currbatch <- currbatch + 1
    blocks[nextsimj:(nextsimj+nhsds-1), c('batch', 'lab', paste0(rep(c('het.b','het.e','het.p'), each = 3), c('','.sd','.cor')) ) := .(currbatch, paste0('all route cor=', cr), 
                                                                                                                                       T, hsds, cr, T, hsds, cr, T, hsds, cr)]
}

CJ.dt = function(X,Y) {
  stopifnot(is.data.table(X),is.data.table(Y))
  k = NULL
  X = X[, c(k=1, .SD)]
  setkey(X, k)
  Y = Y[, c(k=1, .SD)]
  setkey(Y, NULL)
  X[Y, allow.cartesian=TRUE][, k := NULL][]
}
blocks$acute.sc <- NULL
blocksg <- CJ.dt(data.table(country = countries), blocks)
blocksg <- CJ.dt(data.table(acute.sc = acutes), blocksg)
blocksg
blocksg[,c('group','s.epic','s.demog','scale.by.sd','scale.adj','infl.fac','maxN','sample.tmar','psNonPar','each'):= .(country,country, country, T, 1, 200, 10^5, F, F, each.val)]

addParm <- function(x, parmsMat,ii) {
    for(pp in 1:length(parmsMat)) {
        tempP <- as.data.frame(parmsMat)[,pp]
        isch <- !is.numeric(tempP[1])
        parmAdd <- tempP[parmsMat$simNum==ii]
        addStrg <- paste0(" ", names(parmsMat)[pp], "=", "\""[isch], parmAdd, "\""[isch])
        x <- paste0(x, addStrg)
    }
    return(x)
}

blocksg[,simn:=1:nrow(blocksg)]

sink(paste0("HetCounterFactualAcute.txt"))
for(ii in 1:4) { #blocksg[,simn]) {
    cmd <- "R CMD BATCH '--no-restore --no-save --args"
    cmd <- addParm(cmd, blocksg, ii)
    cmd <- paste0(cmd, " ' SimulationStarter.R ", file.path(batchdirnm,'Routs', paste0(nmtmp, sprintf("%06d", ii),'.Rout')), 
                  sep='')
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()

sink("HetCounterFactualAcute.txt")         # create a control file to send to the cluster
if(!file.exists(outdir))      dir.create(outdir) # create directory if necessary
for(aa in acutes)  {                    # loop through acute phase relative hazard
  acdirnm <- file.path(outdir,paste0('Acute', aa)) # set up directory for each acute phase relative hazard
  if(!file.exists(acdirnm))      dir.create(acdirnm) # create directory if necessary
  for(cc in countries) {
    batchdirnm <- file.path(acdirnm, ds.nm[cc]) # setup directory for country
    if(!file.exists(batchdirnm))      dir.create(batchdirnm) # created if necessary
    if(!file.exists(file.path(batchdirnm,'routs')))      dir.create(file.path(batchdirnm, 'routs')) # setup directory to store Rout files
######################################################################
######################################################################           
    ## Set defaults for all parameters for each simulation, simulatin specific-values set later
######################################################################
    ## LEFT OFF HERE
    for(ii in 1:max(sel)) {
      jb <- ii                   # job num
      totn <- totn+1             # total jobs
      cmd <- paste("R CMD BATCH '--args jobnum=", totn, " simj=", ii, " batchdirnm=\"", batchdirnm, "\"", " nc=", nc,
                   " group.ind=", group[ii], " substitute=", substitute, " sub.betas=", sub.betas, " counterf.betas=", counterf.betas,
                   " s.epic=", s.epic[ii],  " s.demog=", s.demog[ii],
                   " s.bmb=", s.bmb[ii], " s.bfb=", s.bfb[ii],
                   " s.bme=", s.bme[ii], " s.bfe=", s.bfe[ii],
                   " s.bmp=", s.bmp[ii], " s.bfp=", s.bfp[ii], 
                   " death=", death[ii],
                   " acute.sc=", aa, " late.sc=", late.sc[ii]," aids.sc=", aids.sc[ii], # acute phase varying throughout loop
                   " bmb.sc=", bmb.sc[ii], " bfb.sc=", bfb.sc[ii],
                   " bme.sc=", bme.sc[ii], " bfe.sc=", bfe.sc[ii],
                   " bmp.sc=", bmp.sc[ii], " bfp.sc=", bfp.sc[ii],
                   " het.b=", het.b[ii], " het.b.sd=", het.b.sd[ii], " het.b.cor=", het.b.cor[ii],
                   " het.e=", het.e[ii], " het.e.sd=", het.e.sd[ii], " het.e.cor=", het.e.cor[ii],
                   " het.p=", het.p[ii], " het.p.sd=", het.p.sd[ii], " het.p.cor=", het.p.cor[ii],                     
                   " het.gen=", het.gen[ii], " het.gen.sd=", het.gen.sd[ii], " het.gen.cor=", het.gen.cor[ii],
                   " het.beh=", het.beh[ii], " het.beh.sd=", het.beh.sd[ii], " het.beh.cor=", het.beh.cor[ii],
                   " hilo=F phihi=.2 phi.m=.2 phi.f=.2 rrhi.m=10 rrhi.f=10", ## default hilo parameters, not used
                   " scale.by.sd=", scale.by.sd[ii], " scale.adj=", scale.adj[ii],
                   " infl.fac=", infl.fac[ii], " maxN=", maxN[ii], " sample.tmar=", sample.tmar[ii],
                   " psNonPar=", psNonPar[ii], " seed=1 tmar=(65*12):(113*12) each=", each[ii],
                   " tint=113*12' SimulationStarter.R ", file.path(batchdirnm, "routs", paste0(ds.nm[group[ii]], ii, ".Rout")), sep='')
  #     if(totn %in% jtd & ii %in% 89:92 & aa==7) { ## for finishing up jobs that didn't get properly submitted (cluster issues sometimes)
 #     if(ii %in% c(1:26,89:92,117:120,145:148) & aa==7) { ## for finishing up jobs that didn't get properly submitted (cluster issues sometimes)
#      if(paste0(group[ii],'-',ii) %in% jtd) {
          num.doing <- num.doing+1
          cat(cmd)               # add command
          cat('\n')              # add new line
 #     }
#}
  }
} 
}
sink()
blocks
totn
save(blocks, file = file.path(outdir,'blocks.Rdata')) # these are country-acute phase specific blocks
####################################################################################################
print(totn)
print(num.doing)
