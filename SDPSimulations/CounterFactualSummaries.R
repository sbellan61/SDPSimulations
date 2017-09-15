####################################################################################################
## Collect, summarize, & visualize results from counterfactual simulations.
####################################################################################################
rm(list=ls())                           # clear workspace
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/DHSProject/SDPSimulations/')
if(grepl('stevenbellan', Sys.info()['login'])) setwd('~/Documents/R Repos/SDPSimulations/SDPSimulations/')
if(grepl('nid', Sys.info()['nodename'])) setwd("/home1/02413/sbellan/SDPSimulations/SDPSimulations")
load("../DHSFitting/data files/dframe.s.Rdata") # country names
source('PlotFunctions.R')                    # load functions to collect & plot results
source('SimulationFunctions.R')                   # load simulation functions
show.pts <- T                           # show observed SDP in real DHS data
## source('CounterFactualSummaries.R')
dir.results <- file.path('results','CounterFactualAll') # results locations
dir.figs <- file.path(dir.results, 'Figures')        # make a directory to store figures
if(!file.exists(dir.figs)) dir.create(dir.figs) # Create directory
load(file.path(dir.results, 'blocksg.Rdata')) ## information on simulations from MK file

resFile <- file.path(dir.results, 'tss.Rdata')

if(file.exists(resFile)) {
    load(resFile)
}else{
    tss <- coll(dir.results, nc = 48, browse=F) ## assemble all results into one DT
    save(tss, file = resFile)
}

tss[,tmp:=paste0(jobnum,'-',yr)] ## Any duplicated jobs X time points?
print(paste(tss[,sum(duplicated(tmp))], 'duplicated jobs'))
tss <- tss[!duplicated(tmp)] ## Remove them
print(paste(tss[,length(unique(jobnum))], 'unique jobs'))
tss$tmp <- NULL
acutes <- tss[,unique(acute.sc)] # which acute phase relative hazards (RHs) were simulated
nac <- length(acutes)                   # how many
countries <- tss[,unique(country)]
countries <- countries[order(countries)]
ncountries <- length(countries)         
jtd <- blocksgTD[!jobnum %in% tss$jobnum, jobnum]
print(paste("didn't do jobs:",paste(head(jtd,50), collapse=','),', ...')) # check to see if any jobs didn't complete
save(jtd, file=file.path(dir.results,'CFJobsToDo.Rdata'))

####################################################################################################
## Figure 2 for the manuscript
####################################################################################################

## matlayout <- t(matrix(c(1:6,10,7:12,10),7,2))

parmsFxn <- function(country=15, s.demog=NA, death=T, acute.sc=5
                   , bmb.sc=1, bme.sc=1, bmp.sc=1
                   , het.b.sd=0, het.b.cor=0, het.e.sd=0, het.e.cor=0, het.p.sd=0, het.p.cor=0
                   , het.beh.sd=0, het.beh.cor=0, het.gen.sd=0, het.gen.cor=0
                   ) {
    if(is.na(s.demog)) s.demog <- country
    parms <- expand.grid(as.list(environment()))
    return(data.table(parms))
}
parmsFxn()
parmsFxn(country=15, death = c(T,F))
nms <- names(parmsFxn())

mains <- c('mortality','pre-couple \ncontact coefficient','extra-couple \ncontact coefficient',
           'intrinsic \ntransmission rate', 'heterogeneity \nin transmission', 'heterogeneity \nwith assortativity')

## choose result contrasts
cc <- 15
sel <-  parmsFxn(country=cc, ) ## base scenario
bases <- tss[sel, .SD, nomatch=0L, on=nms, .SDcols=names(tss)][,unique(jobnum)] ## have duplicates of base simulations based on old code
dups <- bases[-1]
base <- bases[1] ## single base scenario
sel1 <- cbind(parmsFxn(country=cc, death = c(T,F)), cols = c('black','gray'), ltys = c(1,2)) ## mortality contrast
cols <- c('gray','black','brown')
sel2 <- cbind(parmsFxn(country=cc, bmb.sc = c(0,1,10)), cols = cols, ltys = c(2,1,2)) ## pre
sel3 <- cbind(parmsFxn(country=cc, bme.sc = c(0,1,10)), cols = cols, ltys = c(2,1,2)) ## extra
sel4 <- cbind(parmsFxn(country=cc, bmp.sc = c(.1,1,10)), cols = cols, ltys = c(2,1,2)) ## within
cols <- c('black','gray','brown')
sel5 <- cbind(parmsFxn(country=cc, het.gen.sd = c(0,1,2)), cols = cols, ltys = c(1,2,2)) ## gen het
sel6 <- cbind(parmsFxn(country=cc, het.gen.sd = 2, het.gen.cor = c(0,.4,.8)), cols = cols, ltys = c(1,2,2)) ## gen het with assort

## sel5a <- cbind(parmsFxn(country=cc, country=1:16, het.gen.sd = c(0,1,2)), cols = cols, ltys = c(1,2,2)) ## gen het
## sel6a <- cbind(parmsFxn(country=cc, country=1:16, het.gen.sd = 2, het.gen.cor = c(0,.4,.8)), cols = cols, ltys = c(1,2,2)) ## gen het with assort

## xtabs(~country + het.gen.sd, tss[sel5a, .SD, nomatch=0L, on=nms, .SDcols=names(tss)][yr==2013])
## xtabs(~country + het.gen.cor, tss[sel6a, .SD, nomatch=0L, on=nms, .SDcols=names(tss)][yr==2013])

## Make legends
legFxn <- function(ii, leg.cex = .8, cols, ltys
                   ) {
    if(ii==1) legend('bottomleft', c('as fitted', 'no mortality'), lwd = 2, lty = 1, col = cols, bty = 'n', cex = leg.cex)
    if(ii%in%2:3) legend('bottomleft', c('set to 0', 'as fitted', 'scaled X 10'), lwd = 2, lty = ltys,
                         col = cols, bty = 'n', cex = leg.cex)
    if(ii==4) legend(x=1988, y = .5, c('scaled X 1/10', 'as fitted', 'scaled X 10'), lwd = 2, lty = ltys,
                       col = cols, bty = 'n', cex = leg.cex)
    if(ii==5) legend('bottomleft', c('as fitted', 'std dev = 1', 'std dev = 2'), lwd = 2, lty = ltys,
                     col = cols, bty = 'n', cex = leg.cex)
    if(ii==6) legend('bottomleft', c('ro = 0', 'ro = 0.4', 'ro = 0.8'), lwd = 2, lty = ltys,
                     title = 'std dv = 2',
                     col = cols, bty = 'n', cex = leg.cex)
}


cx <- .8
col.pl <- 'black'
pdf(file.path(dir.figs, 'Fig 2 - Counterfactual Summary.pdf'), width = 6.5, h = 5)
matlayout <- t(matrix(c(1:6,13,7:12,13),7,2))
layout(matlayout,w = c(rep(1,6),.8))
par(mar = c(3,1,2,0), oma = c(1,3,0,0), cex.lab = cx, cex.axis = cx, cex.main = cx, fg = col.pl, col.axis = col.pl,
    col.lab = col.pl, col = col.pl, col.main = col.pl)
## par(mfrow = c(2,6), mar = c(3,1,2,.3), oma = c(1,3,0,0), cex.main = .7)
for(ii in 1:6) {
    seltmp <- get(paste0('sel',ii))
    tmp <- tss[seltmp, .SD, nomatch=0L, on=nms, .SDcols=names(tss)][!jobnum %in% dups]
    tmp <- merge(tmp, seltmp) ## add cols, ltys back in
    tmp[, plot(yr, sdp, type = 'n', xlab = '', ylab = '', bty = 'n', xlim=c(1990, 2015), ylim = c(0,1), main =mains[ii], las = 2, yaxt='n')]
    if(ii==1) axis(2, at = seq(0,1,l=5), las = 2) else axis(2, at = seq(0,1,l=5), labels = NA)
    tmp[, lines(yr, sdp, col=cols[1], lty = ltys[1], lwd=2), jobnum]
    points(dframe.s$yr[dframe.s$group==cc], dframe.s$psdc[dframe.s$group==cc], pch = 19, col = 'black', cex = 1.5) 
    seltmp[, legFxn(ii, cols=cols, ltys = ltys)]
}
title(ylab='serodiscordant proportion', outer=T, adj = .5, line = 2) 
####################################################################################################
## Second row
logd <- rep('',6)
logd[2:4] <- 'x'
xlims <- list(c(-.2, 1.2), c(.04, 10.5), c(-.2,3.2), c(-.1, .9))
xlims <- xlims[c(1,2,2,2,3,4)]
xlabs <- c(rep('',4), expression(sigma), expression(rho))
ylab <- ''
for(ii in 1:6) {
    plot(1,1, type = 'n', xlim = xlims[[ii]], ylim = c(0,1), bty = 'n', axes = F, 
         xlab = '', ylab = ylab, main = '', log = logd[ii])
    title(xlab=xlabs[ii], line = 2)
    if(ii==1)     axis(1, at = c(0,1), c('no AIDS \nmortality','as fitted'), las =1, padj = 1)
    if(ii %in% c(2:4))    {
        axis(1, at = c(.1,.2,.5,1,2,5,10), label = c('0.1','0.2','0.5','1','2','5','10'), las = 2)
        if(ii<4) axis(1, at = c(.05), '0', las = 2)
    }
    if(ii==5)             axis(1, at = 0:3)
    if(ii==6)             axis(1, at = seq(0,.8, by = .4))
    if(ii==1) axis(2, at = seq(0,1,l=5), las = 2) else axis(2, at = seq(0,1,l=5), labels = NA)
    if(ii<5) abline(v=1, lwd = 2, col = 'gray')
    if(ii==5) abline(v=0, lwd = 2, col = 'gray')
    for(cc in countries) {
        sel <-  parmsFxn(country=cc) ## base scenario
        bases <- tss[sel, .SD, nomatch=0L, on=nms, .SDcols=names(tss)][,unique(jobnum)] ## have duplicates of base simulations based on old code
        dups <- bases[-1]
        base <- bases[1] ## single base scenario
        ## choose result contrasts
        if(ii==1) sel <- parmsFxn(country=cc, death = c(T,F))
        if(ii==2) sel <- parmsFxn(country=cc, bmb.sc = tss[,unique(bmb.sc)])
        if(ii==3) sel <- parmsFxn(country=cc, bme.sc = tss[,unique(bme.sc)])
        if(ii==4) sel <- parmsFxn(country=cc, bmp.sc = tss[,unique(bmp.sc)])
        if(ii==5) sel <- parmsFxn(country=cc, het.gen.sd = tss[,unique(het.gen.sd)])
        if(ii==6) sel <- parmsFxn(country=cc, het.gen.sd = 2, het.gen.cor = tss[,unique(het.gen.cor)])
        tmp <- tss[sel, .SD, nomatch=0L, on=nms, .SDcols=names(tss)][!jobnum %in% dups]
        ## when duplicate runs were done with same parms (could cause trouble
        ## later if changing the parameters in the sensitivity analysis
        tmp <- unique(tmp[yr==2008,.(death, bmb.sc, bme.sc, bmp.sc, het.gen.sd, het.gen.cor, sdp, yr)]) 
        tmp[, class(het.gen.sd)]
        if(ii==1) tmp[yr==2008, lines(death, sdp, col = dframe$col[cc])]
        if(ii==2) tmp[yr==2008, lines(bmb.sc, sdp, col = dframe$col[cc])]
        if(ii==3) tmp[yr==2008, lines(bme.sc, sdp, col = dframe$col[cc])]
        if(ii==4) tmp[yr==2008, lines(bmp.sc, sdp, col = dframe$col[cc])]
        if(ii==5) tmp[yr==2008][order(het.gen.sd)][, lines(het.gen.sd, sdp, col = dframe$col[cc])]
        if(ii==6) tmp[yr==2008, lines(het.gen.cor, sdp, col = dframe$col[cc])]
    }
}
title(xlab='scalar multiple of fitted parameter used', outer = T, adj = .4, line = -.5)
par(mar=c(3,0,0,0))
plot(0,0, axes = F, bty='n', type = 'n')
legend('bottom', leg = dframe$country, pch = 16, col = dframe$col, cex = .8)
graphics.off()

 
