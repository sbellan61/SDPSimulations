####################################################################################################
## Makes control files for each analysis within which each line giving one R CMD BATCH command line
## to run on a cluster.
####################################################################################################
rm(list=ls())                                  # clear workspace
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/DHSProject/SDPSimulations/')
if(grepl('stevenbellan', Sys.info()['login'])) setwd('~/Documents/R Repos/SDPSimulations/SDPSimulations/')
source("SimulationFunctions.R")                   # load simulation functions from script
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') #  transmission coefficient names, for convenience
nc <- 24                                       # core per simulation
## source('AssortHiLoMK.R')

####################################################################################################
####################################################################################################
##   Assortativity COUNTERFACTUAL ANALYSIS
####################################################################################################
## Using an HIV transmission model fit to Demographic and Health Surveys in 18 sub-Saharan African
## countries, we used counter-factual simulations to examine how serodiscordant proportions are
## affected by assortativity between patners
####################################################################################################
## countries <- 1:length(ds.nm)
countries <- which(ds.nm=='Zambia')
each <- 200 ##  equates to ~100,000 couples
blocks <- expand.grid(acute.sc = 7,
                      late.sc = 1, aids.sc = 1, death = T,
                      hilo = T,
                      rel.assort = c(1:3), ## multiplier of phi^2 to tune assortativity
                      phi = c(.01, .05, .1, .3), ## proportion of population that's high risk
                      rrhi.m = c(1, 5, 10, 50), ## relative risk high vs low
                      bmb.sc = 1, ## bfb.sc = 1,
                      bme.sc = 1, ## bfe.sc = 1,
                      bmp.sc = c(1,2,5), 
                      country = countries)
blocks$bfb.sc <- blocks$bmb.sc
blocks$bfe.sc <- blocks$bme.sc
blocks$bfp.sc <- blocks$bmp.sc
blocks$phihi <- blocks$phi^2*blocks$rel.assort
blocks$phi.m <- blocks$phi - blocks$phi^2*blocks$rel.assort
blocks$phi.f <- blocks$phi.m
blocks$rrhi.f <- blocks$rrhi.m
blocks$jobnum <- 1:nrow(blocks)
blocksg <- as.data.table(blocks)

blocks

with(blocks, phihi + phi.m)

out.dir <- file.path('results','AssortHetHeatMap')
if(!file.exists(out.dir))      dir.create(out.dir) # create directory if necessary
if(!file.exists(file.path(out.dir,'Rdatas')))      dir.create(file.path(out.dir,'Rdatas')) # create directory if necessary
if(!file.exists(file.path(out.dir,'Routs')))      dir.create(file.path(out.dir,'Routs')) # create directory if necessary

blocksg[,c('group','s.epic','s.demog','scale.by.sd','scale.adj','infl.fac','maxN','sample.tmar','psNonPar','each'):= .(country,country, country, T, 1, 200, 10^5, F, F, 200)]
blocksg[,jobnum:=1:nrow(blocksg)]
blocksg[,c('seed','out.dir','sim.nm','doSubs', 'sub.betas'):=.(1,out.dir, 'AHL', F,F)]
blocksg[,c('tmar','tint'):=.('tmar=(65*12):(113*12)',113*12)]
blocksg[,c('acute.sc'):=.(5)]

blocksgTD <- blocksg[country==15] ##

sink("AssortHiLo.txt") ## create a control file to send to the cluster
for(ii in blocksgTD[,jobnum]) { #blocksgTD[,jobnum]) {
    cmd <- "R CMD BATCH '--no-restore --no-save --args"
    cmd <- addParm(cmd, blocksgTD, ii) ## remove lab since it has spaces & isn't used in psrun
    cmd <- paste0(cmd, " ' SimulationStarter.R ", file.path(out.dir,'Routs', paste0('AHLsim', sprintf("%06d", ii),'.Rout')), 
                  sep='')
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()
save(blocksg, file = file.path(out.dir,'blocksg.Rdata')) # these are country-acute phase specific blocks
print(nrow(blocksg))
print(nrow(blocksgTD))
