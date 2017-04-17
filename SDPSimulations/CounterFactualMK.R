####################################################################################################
## Makes control files for each analysis within which each line giving one R CMD BATCH command line
## to run on a cluster.
####################################################################################################

## source('CounterFactualMK.R')
#rm(list=ls())                                  # clear workspace
require(data.table)
source("SimulationFunctions.R")                   # load simulation functions from script
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

out.dir <- file.path('results','CounterFactual')

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
blocksg[,jobnum:=1:nrow(blocksg)]
blocksg[,c('seed','out.dir','sim.nm','doSubs'):=.(1,out.dir, 'CF', F)]
blocksg[,c('tmar','tint'):=.('tmar=(65*12):(113*12)',113*12)]

blocksgTD <- blocksg[country==15]

if(!file.exists(out.dir))      dir.create(out.dir) # create directory if necessary
if(!file.exists(file.path(out.dir,'Rdatas')))      dir.create(file.path(out.dir,'Rdatas')) # create directory if necessary
if(!file.exists(file.path(out.dir,'Routs')))      dir.create(file.path(out.dir,'Routs')) # create directory if necessary
sink("HetCounterFactualAcute.txt")         # create a control file to send to the cluster
for(ii in blocksgTD[,jobnum]) { #blocksgTD[,jobnum]) {
    cmd <- "R CMD BATCH '--no-restore --no-save --args"
    cmd <- addParm(cmd, blocksgTD[, !'lab'], ii) ## remove lab since it has spaces & isn't used in psrun
    cmd <- paste0(cmd, " ' SimulationStarter.R ", file.path(out.dir,'Routs', paste0('CFsim', sprintf("%06d", ii),'.Rout')), 
                  sep='')
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()
save(blocksg, file = file.path(out.dir,'blocksg.Rdata')) # these are country-acute phase specific blocks
print(nrow(blocksg))
print(nrow(blocksgTD))
