####################################################################################################
## Collect, summarize, & visualize results from simulations of one
## country with varying behavior, assortativity & heterogeneity.
####################################################################################################
rm(list=ls())                           # clear workspace
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/DHSProject/SDPSimulations/')
source('PlotFunctions.R')                    # load functions to collect & plot results
source('SimulationFunctions.R')                   # load simulation functions
library(abind)                          # array binding
load('../DHSFitting/data files/ds.nm.all.Rdata')        # country names
load('../DHSFitting/data files/dframe.s.Rdata')        # SDP by survey
load('../DHSFitting/data files/dframe.Rdata') # country summary data (peak prevalence, country-prevalence, etc...)
load('../DHSFitting/data files/draw.Rdata') # raw country summary data (peak prevalence, country-prevalence, etc...)
load('../DHSFitting/data files/dframe.s.Rdata')# country summary data (peak prevalence, country-prevalence, etc...) by DHS survey
load('../DHSFitting/data files/draw.s.Rdata')# country summary data (peak prevalence, country-prevalence, etc...) by DHS survey (unfiltered data)
do.again <- T                           # collect results again (otherwise load cfs.Rdata)
show.pts <- T                           # show observed SDP in real DHS data
## source('CountryAssortSimSummaries.R')

dir.results <- file.path('results','AssortBehaviorCountrySim')
dir.figs <- file.path(dir.results, 'Figures')        # make a directory to store figures
if(!file.exists(dir.figs)) dir.create(dir.figs)      # create it
## load 'blocks' which gives info on all simulations within each country-acute group
load(file.path(dir.results, 'blocks.Rdata'))
## Load all results files in the results directory (all Rdata files except blocks & cfs)
fs <- list.files(pattern = '.Rdata', path = file.path(dir.results), recursive = T, full.names = T)
fs <- fs[!grepl('blocks',fs) & !grepl('cfs',fs) &  !grepl('JobsToDo',fs)]
#fs <- fs[grepl('Acute7',fs)]
print(length(fs))

## coll: collects all results into data frame of input parameters, and array of time series
if(!file.exists(file.path(dir.results, 'cfs.Rdata')) | do.again) {
  cfs <- coll(fs, nc = 12, give.ev = F, lbrowse=F, trace=F, browse=F)
  cfs$cframe <- cfs$cframe[order(cfs$cframe$job),] # order by jobs
  cfs$t.arr <- cfs$t.arr[,,order(cfs$cframe$job)]          # ditto
  attach(cfs) # for convenience, be careful later!
  save(cfs, file = file.path(dir.results, 'cfs.Rdata'))
}else{
  load(file.path(dir.results, 'cfs.Rdata'))
  attach(cfs)
}
dim(cframe) ## dimensions & check for duplicate runs
print(paste(sum(duplicated(cframe$job)), 'duplicate jobs'))
cframe <- cframe[!duplicated(cframe$job),]
t.arr <-  t.arr[,,!duplicated(cframe$job)]
acutes <- unique(cframe$acute.sc)       # which acute phase relative hazards (RHs) were simulated
nac <- 8 #length(acutes)                   # how many
col.pl <- 'black'                       # base plot color
mend <- max(blocks$end)                 # last job of each block
countries <- unique(cframe$group.ind)   # countries simulated
countries <- countries[order(countries)]
ngroup <- length(countries)             # how many
jtd.all <- with(blocks, jobnum[group==12])
jtd <- jtd.all[!jtd.all %in% cframe$job]
print(paste("didn't do", length(jtd), "jobs")) # check to see if any jobs didn't complete
save(jtd, file=file.path(dir.results, 'JobsToDo.Rdata'))

temp <- cframe[,c('het.beh.sd','m.het.beh.var','f.het.beh.var', 'het.beh.cor','het.beh.cov.emp')]
temp$het.beh.cor <- temp$het.beh.cor*temp$het.beh.sd
unique(temp)

het.sd <- 2
het.cor <- .7
b.corSig <- matrix(c(het.sd^2,rep(het.cor*het.sd^2,2),het.sd^2),2,2)
hets <- rmvnorm(10^4, mean = rep(0,2), sigma = b.corSig)
b.corSig
cov(hets)

####################################################################################################
## Need to fit these simulations with a homogenous model

 
####################################################################################################
## Plot SDP vs prev
prev.max <- .35
sdp.mins <- c(0,.35)
sdp.maxs <- c(1,.9)
bcors <- round(unique(blocks$het.beh.cor)^2,2)
bsds <- unique(blocks$het.beh.sd)
cframe$het.beh.cor2 <- round(cframe$het.beh.cor^2,2)
cols <- colorRampPalette(c('yellow','red'))(21)
for(ii in 1:2) {
    pdf(file.path(dir.figs, paste0('SDP vs Prev', c('full','zoom')[ii],'.pdf')), w = 6.5, h = 6)
    layout(cbind(t(matrix(1:16, 4,4)), rep(17,4)))
    par(mar = c(2,2,.5,.5), oma = c(6,7,0,0), 'ps' = 12)
    for(cr in rev(bcors)) {
        for(sd in bsds) {
            if(sd==0) crr <- 0 else crr <- cr
                                        #sd <- 1; crr = .75
            sel <- which(cframe$het.beh.sd==sd & cframe$het.beh.cor2 == crr)
            c2 <- cframe[sel,c('bme.sc','sdp08','prev08','het.beh.sd','het.beh.cor')]
            plot(1,1, type = 'n', las = 1, bty = 'n', ylim = c(sdp.mins[ii],sdp.maxs[ii]), xlim = c(0,prev.max), xlab = '', ylab = '', axes = F)
            with(dframe, points(dhsprev, psdc, pch = 16, col = gray(.6)))
            with(cframe, points(prev08[sel], sdp08[sel], col = cols, pch = 16))
            if(sd==bsds[1]) {
                axis(2, pretty(c(sdp.mins[ii],sdp.maxs[ii]),5), las = 2)
                mtext(cr, 2, adj = .5, line = 7, las = 2)
                if(cr==bcors[4]) mtext(expression(rho['contact']), 2, at = .9, line = 5, las = 2)
            }else{
                axis(2, pretty(c(sdp.mins[ii],sdp.maxs[ii]),5), lab = NA)
            }
            if(cr==bcors[1]) {
                axis(1, pretty(c(0,prev.max),5))
                if(sd==0) mtext(expression(sigma['contact']), 1, at = -.2, line = 7)
                mtext(sd, 1, adj = .5, line = 7)
            }else{
                axis(1, pretty(c(0,prev.max),5), lab = NA)
            }
        }
    }
    mtext('SDP', 2 , 2, T)
    mtext('prevalence', 1 , 2, T)
    par(mar=rep(0,4))
    plot(0,0, type = 'n', bty = 'n', axes = F)
    legend('center', leg = signif(cframe$bmb.sc[sel],2), col = cols, pch = 16, title = 'scalar multiplier \nof fitted\n contact coefficients', bty ='n')
    legend('top', leg = 'DHS country data', col = gray(.6), pch = 16, bty = 'n')
    graphics.off()
}
xtabs(~het.beh.sd + het.beh.cor, cframe)
