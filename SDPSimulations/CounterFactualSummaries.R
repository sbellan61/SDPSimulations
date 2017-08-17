####################################################################################################
## Collect, summarize, & visualize results from counterfactual simulations.
####################################################################################################
rm(list=ls())                           # clear workspace
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/DHSProject/SDPSimulations/')
if(grepl('stevenbellan', Sys.info()['login'])) setwd('~/Documents/R Repos/SDPSimulations/SDPSimulations/')
load("../DHSFitting/data files/dframe.s.Rdata") # country names
source('PlotFunctions.R')                    # load functions to collect & plot results
source('SimulationFunctions.R')                   # load simulation functions
show.pts <- T                           # show observed SDP in real DHS data
## source('CounterFactualSummaries.R')
dir.results <- file.path('results','CounterFactual') # results locations
dir.figs <- file.path(dir.results, 'Figures')        # make a directory to store figures
if(!file.exists(dir.figs)) dir.create(dir.figs) # Create directory
load(file.path(dir.results, 'blocksg.Rdata')) ## information on simulations from MK file
tss <- coll(dir.results, nc = 48, browse=F) ## assemble all results into one DT
save(tss, file = file.path(dir.results, 'tss.Rdata'))

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
jtd <- blocksg[!jobnum %in% tss$jobnum, jobnum]
print(paste("didn't do jobs:",paste(head(jtd,50), collapse=','),', ...')) # check to see if any jobs didn't complete
save(jtd, file=file.path(dir.results,'CFJobsToDo.Rdata'))

####################################################################################################
## Figure 2 for the manuscript
####################################################################################################


## matlayout <- t(matrix(c(1:6,10,7:12,10),7,2))
## layout(matlayout,w = c(rep(1,6),.85))
## par(mar = c(3,1,2,0), oma = c(1,3,0,0), cex.lab = cx, cex.axis = cx, cex.main = cx, fg = col.pl, col.axis = col.pl,
##     col.lab = col.pl, col = col.pl, col.main = col.pl)

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

sel <-  parmsFxn()
bases <- tss[sel, .SD, nomatch=0L, on=nms, .SDcols=names(tss)][,unique(jobnum)] ## have duplicates of base simulations based on old code
dups <- bases[-1]
base <- bases[1]
sel1 <- cbind(parmsFxn(death = c(T,F)), cols = c('black','gray'), ltys = c(1,2))
cols <- c('gray','black','brown')
sel2 <- cbind(parmsFxn(bmb.sc = c(0,1,10)), cols = cols, ltys = c(2,1,2))
sel3 <- cbind(parmsFxn(bme.sc = c(0,1,10)), cols = cols, ltys = c(2,1,2))
sel4 <- cbind(parmsFxn(bmp.sc = c(.1,1,10)), cols = cols, ltys = c(2,1,2))
cols <- c('black','gray','brown')
sel5 <- cbind(parmsFxn(het.gen.sd = c(0,1,2)), cols = cols, ltys = c(1,2,2))
sel6 <- cbind(parmsFxn(het.gen.sd = 2, het.gen.cor = c(0,.4,.8)), cols = cols, ltys = c(1,2,2))


legFxn <- function(ii, leg.cex = .8, cols, ltys
                   ) {
    if(ii==1) legend('bottomleft', c('as fitted', 'no AIDS mortality'), lwd = 2, lty = 1, col = cols, bty = 'n', cex = leg.cex)
    if(ii%in%2:3) legend('bottomleft', c('set to 0', 'as fitted', 'scaled X 10'), lwd = 2, lty = ltys,
                         col = cols, bty = 'n', cex = leg.cex)
    if(ii==4) legend('right', c('scaled X 1/10', 'as fitted', 'scaled X 10'), lwd = 2, lty = ltys,
                       col = cols, bty = 'n', cex = leg.cex)
    if(ii==5) legend('bottomleft', c('as fitted', 'std dev = 1', 'std dev = 2'), lwd = 2, lty = ltys,
                     col = cols, bty = 'n', cex = leg.cex)
    if(ii==6) legend('bottomleft', c('ro = 0', 'ro = 0.4', 'ro = 0.8'), lwd = 2, lty = ltys,
                     title = 'std dv = 2',
                     col = cols, bty = 'n', cex = leg.cex)
}

cc <- 15
pdf(file.path(dir.figs, 'Fig 2 - Counterfactual Summary.pdf'))
par(mfrow = c(2,3))
for(ii in 1:6) {
    tmp[, plot(yr, sdp, type = 'n', xlab = '', ylab = '', bty = 'n', xlim=c(1990, 2015), ylim = c(0,1), main =mains[ii], las = 2)]
    seltmp <- get(paste0('sel',ii))
    tmp <- tss[seltmp, .SD, nomatch=0L, on=nms, .SDcols=names(tss)][!jobnum %in% dups]
    tmp <- merge(tmp, seltmp) ## add cols, ltys back in
    tmp[, lines(yr, sdp, col=cols[1], lty = ltys[1], lwd=2), jobnum]
    points(dframe.s$yr[dframe.s$group==cc], dframe.s$psdc[dframe.s$group==cc], pch = 19, col = 'black', cex = 1.1) 
    seltmp[, legFxn(ii, cols=cols, ltys = ltys)]
}
graphics.off()
 
####################################################################################################

####################################################################################################
### WORKING DOWN HERE
####################################################################################################

        ## Extract SDP's from 2008 for all scenarios to plot in 2nd row of figure
pdf(file.path(dir.figs, 'test.pdf'))
par(mar = c(3,1,.5,0))
cols <- rainbow(length(ds.nm))
ylab <- '' #ifelse(rr==1,'SDP','')
xmax <- c(1,10,10,10,2)
for(rr in 1:4) {
    yaxt <- ifelse(rr==1, T, F)
    if(rr %in% 2:4) {
        logd <- 'x'
        xlim <- c(.04, 10.5)
    }else{
        logd <- ''
        if(rr==1)  xlim <- c(-.2,1.2)     else    xlim <- c(-.2,3.2)
    }   
    plot(1,1, type = 'n', xlim = xlim, ylim = c(0,1), bty = 'n', axes = F,
         xlab = '', ylab = ylab, main = '', log = logd)

yrind <- which(t.arr[,2,1]==2008)    

ac <- 5
cc <- 15
logd <- rep('',6)
logd[2:4] <- 'x'
xlims <- list(c(-.2, 1.2), c(.04, 10.5), c(-.2,3.2))
xlims <- xlims[c(1,2,2,2,3,1)]
ylab <- ''

ii=1

    countryCols <- 'black' #rainbow(length(countries))
    yaxt <- ifelse(ii==1, T, F)
    pdf(file.path(dir.figs, 'test.pdf'))
    for(ii in 1:6) {
        plot(1,1, type = 'n', xlim = xlims[[ii]], ylim = c(0,1), bty = 'n', axes = F,
             xlab = '', ylab = ylab, main = '', log = logd[ii])
        ## x-axis
        if(ii==1)     axis(1, at = c(0,1), c('no AIDS \nmortality','as fitted'), las =1, padj = 1)
        if(ii %in% c(2:4))    {
            axis(1, at = c(.1,.2,.5,1,2,5,10), label = c('0.1','0.2','0.5','1','2','5','10'), las = 2)
            if(ii<4) axis(1, at = c(.05), '0', las = 2)
        }
        if(ii==5)             axis(1, at = 0:3)
        if(ii ==5) {
            grcol <- rgb(t(col2rgb(gray(.6), alpha = .5)),max=255)
            segments(0,0,0,1, col = grcol, lwd = 3)
        }else{
            grcol <- rgb(t(col2rgb(gray(.6), alpha = .5)),max=255)
            segments(1,0,1,1, col = grcol, lwd = 3)
        }
        ## y-axis
        if(yaxt) axis(2, at = seq(0,1,l=5), las = 2) else axis(2, at = seq(0,1,l=5), labels = NA)
        for(cc in countries)   {        ## loop thru countries
            inds <- indsFxn(cc=cc, ac=ac)
            for(jj in 1:length(inds[[ii]])) {
                load(fls[jj])
                tsstmp <- data.table(output$tss)
                tsstmp[, sdp:=(mm+ff)/inf.alive]
                tsstmp <- cbind(tsstmp, data.table(t(output$pars)), simNumplot = jj)
                if(jj==1) {
                    tss <- tsstmp
                }else{
                    tss <- rbind(tss,tsstmp)
                }
            }
            if(ii==1) tss[yr==2008, lines(c(1,0), sdp, col = 'black')] #countryCols[cc])]
            if(ii==2) tss[yr==2008, lines(, sdp, col = 'black')] #countryCols[cc])]
        }
    }
    graphics.off()

        tss[yr==2008, sdp]
xlims[[ii]]

                    if(rr==1) { 


            sel <- which(cframe$group.ind==cc & cframe$simj %in% 1:2 & cframe$acute.sc==ac)
            
            lines(c(1,0), (t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel], col = cols[cc], type = 'l', pch = 19, lty = 1)
        }
        if(rr==2) {
            sel <- which(cframe$group.ind==cc & cframe$simj %in% blocksg$start[bltd[rr]]:blocksg$end[bltd[rr]] & cframe$acute.sc==ac)
            xs <- cframe$bmb.sc[sel]
            xs[cframe$bmb.sc[sel]==0] <- .05 # since on a log scale
            lines(xs, (t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel],
                  col = cols[cc], type = 'l', pch = 19, lty = 1)
        }
        if(rr==3) {
            sel <- which(cframe$group.ind==cc & cframe$simj %in% blocksg$start[bltd[rr]]:blocksg$end[bltd[rr]] & cframe$acute.sc==ac)
            xs <- cframe$bme.sc[sel]
            xs[cframe$bme.sc[sel]==0] <- .05 # since on a log scale
            lines(xs, (t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel],
                  col = cols[cc], type = 'l', pch = 19, lty = 1)
        }
        if(rr %in% 4) {
            sel <- which(cframe$group.ind==cc & cframe$simj %in% (blocksg$start[bltd[rr]]+1):blocksg$end[bltd[rr]] & cframe$acute.sc==ac) 
            lines(cframe$bmp.sc[sel], (t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel], col = cols[cc], type = 'l', pch = 19, lty = 1)
        }
        if(rr == 5) {
            sel <- which(cframe$group.ind==cc & cframe$simj %in% c(1,90:92) & cframe$acute.sc==ac)
            lines(cframe$het.gen.sd[sel], (t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel],
                  col = cols[cc], type = 'l', pch = 19, lty = 1)
        } 
    }
}
mtext('serodiscordant proportion (SDP)', side = 2, line = 2, adj = .5, cex = leg.cex, outer = T)
mtext('scalar multiple of fitted parameter used', side = 1, line = -.3, adj = .35, cex = leg.cex, outer = T)
mtext('standard deviation \nof risk distribution', side = 1, line = 3, adj = .5, cex = leg.cex, outer = F)      
plot.new()
par(mar=c(0,0,0,0))
sel <- which(cframe$simj==1 & cframe$acute.sc==ac)
ord <- rev(order((t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel]))
save(ord, file=file.path(dir.results,'sdp ord.Rdata'))
par(xpd=NA)
legend('bottomleft', ds.nm[ord], col=rainbow(length(sel))[ord], lwd = 1, title = 'top to bottom', cex = .8, inset = -.1)
if(one.file) mtext(ds.nm[cc1], adj = 0.5, line = -2, side = 3, outer = F, cex = .75)
if(!one.file) dev.off()
}
if(one.file)    dev.off()
}
graphics.off()
####################################################################################################
####################################################################################################
####################################################################################################


 
