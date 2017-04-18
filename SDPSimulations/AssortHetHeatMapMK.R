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
## source('AssortHetHeatMapMK.R')

####################################################################################################
## Simulate Tanzania but with assortativity/heterogeneity & varied amounts of pre-couple extra-couple
####################################################################################################
## cc <- which(ds.nm=='Tanzania') 
countries <- 1:length(ds.nm)
each.val <- 200 ##  equates to ~100,000 couples

blocks.beh <- expand.grid(country = countries, het.beh=T, het.beh.sd=seq(0,3, by = .1), het.beh.cor=seq(0,1, by = .05), bmb.sc = 1)
blocks.beh$bfb.sc <- blocks.beh$bme.sc <- blocks.beh$bfe.sc <- blocks.beh$bmb.sc ## all contact coefficients scaled same way
blocks.beh$bmp.sc <- blocks.beh$bfp.sc <- 1
blocks.beh$het.gen <- F; blocks.beh$het.gen.sd <- blocks.beh$het.gen.cor <- 0
blocks.beh <- blocks.beh[with(blocks.beh, order(country, het.beh.sd, het.beh.cor, bmb.sc)),]

blocks.gen <- expand.grid(country = countries, het.gen=T, het.gen.sd=seq(0,3, by = .1), het.gen.cor=seq(0,1, by = .05), bmb.sc = 1)
blocks.gen$bfb.sc <- blocks.gen$bme.sc <- blocks.gen$bfe.sc <- blocks.gen$bmb.sc ## all contact coefficients scaled same way
blocks.gen$bmp.sc <- blocks.gen$bfp.sc <- 1
blocks.gen$het.beh <- F; blocks.gen$het.beh.sd <- blocks.gen$het.beh.cor <- 0
blocks.gen <- blocks.gen[with(blocks.gen, order(country, het.beh.sd, het.beh.cor, bmb.sc)),]

blocks <- rbind(blocks.beh, blocks.gen)
blocks$jobnum <- 1:nrow(blocks)
blocksg <- data.table(blocks)

nn <- nrow(blocks)
head(blocks,50)

out.dir <- file.path('results','AssortHetHeatMap')
blocksg[,c('group','s.epic','s.demog','scale.by.sd','scale.adj','infl.fac','maxN','sample.tmar','psNonPar','each'):= .(country,country, country, T, 1, 200, 10^5, F, F, each.val)]
blocksg[,jobnum:=1:nrow(blocksg)]
blocksg[,c('seed','out.dir','sim.nm','doSubs'):=.(1,out.dir, 'AHH', F)]
blocksg[,c('tmar','tint'):=.('tmar=(65*12):(113*12)',113*12)]
blocksg[,c('acute.sc'):=.(5)]

blocksgTD <- blocksg[country==15] ##

if(!file.exists(out.dir))      dir.create(out.dir) # create directory if necessary
if(!file.exists(file.path(out.dir,'Rdatas')))      dir.create(file.path(out.dir,'Rdatas')) # create directory if necessary
if(!file.exists(file.path(out.dir,'Routs')))      dir.create(file.path(out.dir,'Routs')) # create directory if necessary
sink("AssortHetHeatMap.txt")         # create a control file to send to the cluster
for(ii in blocksgTD[,jobnum]) { #blocksgTD[,jobnum]) {
    cmd <- "R CMD BATCH '--no-restore --no-save --args"
    cmd <- addParm(cmd, blocksgTD, ii) ## remove lab since it has spaces & isn't used in psrun
    cmd <- paste0(cmd, " ' SimulationStarter.R ", file.path(out.dir,'Routs', paste0('AHHsim', sprintf("%06d", ii),'.Rout')), 
                  sep='')
    cat(cmd)               # add command
    cat('\n')              # add new line
}
sink()
save(blocksg, file = file.path(out.dir,'blocksg.Rdata')) # these are country-acute phase specific blocks
print(nrow(blocksg))
print(nrow(blocksgTD))
