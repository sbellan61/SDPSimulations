####################################################################################################
## Collect, summarize, & visualize results from simulations of one
## country with varying assortativity & heterogeneity.
####################################################################################################
library(plotrix)
rm(list=ls())                           # clear workspace
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/DHSProject/SDPSimulations/')
source('PlotFunctions.R')                    # load functions to collect & plot results
source('SimulationFunctions.R')                   # load simulation functions
library(abind)                          # array binding
load('data files/ds.nm.all.Rdata')        # country names
load('data files/dframe.s.Rdata')        # SDP by survey
do.again <- T                           # collect results again (otherwise load cfs.Rdata)
show.pts <- T                           # show observed SDP in real DHS data
## source('AssortHetHeatMapSummary.R')

dir.results <- file.path('results','AssortHetHeatMap')
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
jtd.all <- with(blocks, jobnum[group %in% 11:12])
jtd <- jtd.all[!jtd.all %in% cframe$job]
print(paste("didn't do jobs:",length(jtd), collapse=',')) # check to see if any jobs didn't complete
save(jtd, file=file.path(dir.results, 'JobsToDo.Rdata'))

####################################################################################################
## 
par(mar=c(4,4,.3,0))

xs <- unique(cframe$het.beh.sd)
nx <- length(xs)
ys <- unique(cframe$het.beh.cor)
ny <- length(ys)

## Behavior het
out.beh <- with(cframe[cframe$het.beh==T,], {
    prev <- sdp <- matrix(NA, nx, ny)
    for(ii in 1:nx) {
        for(jj in 1:ny) {
            sel <- which(het.beh.sd==xs[ii] & het.beh.cor==ys[jj])
            if(length(sel)>0) {
                sdp[ii,jj] <- sdp08[sel]
                prev[ii,jj] <- prev08[sel]
            }
        }}
    return(list(sdp,prev))
})
sdp.beh <- out.beh[[1]]
prev.beh <- out.beh[[2]]
## Gen het
out.gen <- with(cframe[cframe$het.gen==T,], {
    prev <- sdp <- matrix(NA, nx, ny)
    for(ii in 1:nx) {
        for(jj in 1:ny) {
            sel <- which(het.gen.sd==xs[ii] & het.gen.cor==ys[jj])
            if(length(sel)>0) {
                sdp[ii,jj] <- sdp08[sel]
                prev[ii,jj] <- prev08[sel]
            }
        }}
    return(list(sdp,prev))
})
sdp.gen <- out.gen[[1]]
prev.gen <- out.gen[[2]]

for(show.prev in c(T,F)) {
    ht <- ifelse(show.prev, 5, 3)
    for(ii in 1:2) {
        if(ii==1)       pdf(file.path(dir.figs, paste0('HeatPlot', ' prev'[show.prev], '.pdf')), w = 6.5, h = ht)
        if(ii==2)       jpeg(file.path(dir.figs, paste0('HeatPlot', ' prev'[show.prev], '.jpg')), w = 6.5, h = ht, unit ='in', res = 300)
        par('ps'=12, mar = c(4,4,1,2))
        if(show.prev)       layout(t(matrix(c(1:6),3,2)), w = c(1,1,.3), h = c(1,1))
        if(!show.prev)       layout(matrix(c(1:3),1,3), w = c(1,1,.3))
        ## ##################################################
        ## SDP
        ## ##################################################
        levels <- seq(.5,1, b = .01)
        cols <- colorRampPalette(c('yellow','red'))(length(levels)-1)
        ## Contact rates
        image(xs, ys, sdp.beh,  xlim = c(-.1,3), ylim = c(-.1,1), mgp = c(3,1,0), bty='n', las = 1,
              col = cols, breaks = levels, axes = T, xlab = expression(sigma['contact']), ylab = expression(rho['contact']))
        ## Transmission Rates
        image(xs, ys, sdp.gen,  xlim = c(-.1,3), ylim = c(-.1,1), mgp = c(3,1,0), bty='n', las = 1,
              col = cols, breaks = levels, axes = T, xlab = expression(sigma['transmission']), ylab = expression(rho['transmission']))
        ## Palette legend
        par(mar=rep(0,4))
        plot(0,0,type="n",axes=F, xlim = c(-.1,.2), ylim = c(-.1,.9), xlab = '', ylab = '')
        sel <- which(1:length(levels) %% 5 == 1)
        color.legend(.09,.1,.15,.8, levels[sel], rect.col = cols[sel-1], gradient = "y", cex = 1)
        text(.025, .8, 'SDP', pos = 3, cex = 1)
        if(show.prev) {
            ## ##################################################
            ## Prevalence
            ## ##################################################
            par('ps'=12, mar = c(4,4,1,2))
            levels <- seq(0,.1, b = .001)
            cols <- colorRampPalette(c('yellow','red'))(length(levels)-1)
            ## Contact rates
            image(xs, ys, prev.beh,  xlim = c(-.1,3), ylim = c(-.1,1), mgp = c(3,1,0), bty='n', las = 1,
                  col = cols, breaks = levels, axes = T, xlab = expression(sigma['contact']), ylab = expression(rho['contact']))
            ## Transmission Rates
            image(xs, ys, prev.gen,  xlim = c(-.1,3), ylim = c(-.1,1), mgp = c(3,1,0), bty='n', las = 1,
                  col = cols, breaks = levels, axes = T, xlab = expression(sigma['transmission']), ylab = expression(rho['transmission']))
            ## Palette legend
            par(mar=rep(0,4))
            plot(0,0,type="n",axes=F, xlim = c(-.1,.2), ylim = c(-.1,.9), xlab = '', ylab = '')
            sel <-  which(1:length(levels) %% 10 == 1)
            color.legend(.09,.1,.15,.8, levels[sel], rect.col = cols[sel-1], gradient = "y", cex = 1)
            text(.025, .8, 'prevalence', pos = 3, cex = 1)
        }
        graphics.off()
    }
}
