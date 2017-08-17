####################################################################################################
## Collect, summarize, & visualize results from counterfactual simulations.
####################################################################################################
rm(list=ls())                           # clear workspace
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/DHSProject/SDPSimulations/')
if(grepl('stevenbellan', Sys.info()['login'])) setwd('~/Documents/R Repos/SDPSimulations/SDPSimulations/')
load("../DHSFitting/data files/dframe.s.Rdata") # country names
source('PlotFunctions.R')                    # load functions to collect & plot results
source('SimulationFunctions.R')                   # load simulation functions
do.again <- F                           # collect results again (otherwise load cfs.Rdata)
show.pts <- T                           # show observed SDP in real DHS data
## source('CounterFactualSummaries.R')
dir.results <- file.path('results','CounterFactual') # results locations
dir.figs <- file.path(dir.results, 'Figures')        # make a directory to store figures

if(!file.exists(dir.figs)) dir.create(dir.figs)      # create it
## load 'blocks' which gives info on all simulations within each country-acute group
load(file.path(dir.results, 'blocksg.Rdata'))
## Load all results files in the results directory (all Rdata files except blocks & cfs)
fs <- list.files(path = file.path(dir.results, 'Rdatas'), recursive = T, full.names = T)
print(length(fs))

## coll: collects all results into data frame of input parameters, and array of time series
if(!file.exists(file.path(dir.results, 'cfs.Rdata')) | do.again) {
    tss <- coll(dir.results, nc = 48, browse=F)
    save(tss, file = file.path(dir.results, 'tss.Rdata'))
}else{
    load(file.path(dir.results, 'tss.Rdata'))
}

tss[,length(unique(jobnum))]
dim(tss) ## dimensions & check for duplicate runs
print(paste(tss[, sum(duplicated(jobnum))]), 'duplicated jobs')
tss <- tss[!duplicated(jobnum)]
acutes <- tss[,unique(acute.sc)] # which acute phase relative hazards (RHs) were simulated
nac <- length(acutes)                   # how many
countries <- tss[,unique(country)]
countries <- countries[order(countries)]
ncountries <- length(countries)         
jtd <- blocksg[!jobnum %in% tss$jobnum, jobnum]
print(paste("didn't do jobs:",paste(head(jtd,50), collapse=','),', ...')) # check to see if any jobs didn't complete
save(jtd, file=file.path(dir.results,'CFJobsToDo.Rdata'))

blocks <- blocksg[jobnum %in% cframe$job, ]
ac.to.do <- c(1) ## 1,7,25,50)
nac <- length(ac.to.do)

####################################################################################################
## Figure 2 for the manuscript
####################################################################################################

## matlayout <- t(matrix(c(1:6,10,7:12,10),7,2))
## layout(matlayout,w = c(rep(1,6),.85))
## par(mar = c(3,1,2,0), oma = c(1,3,0,0), cex.lab = cx, cex.axis = cx, cex.main = cx, fg = col.pl, col.axis = col.pl,
##     col.lab = col.pl, col = col.pl, col.main = col.pl)

indsFxn <- function(cc, ac) {
inds <- list()
length(inds) <- 6
    inds[[1]] <- c(blocks[group==cc & acute.sc==ac & lab=="as fitted", jobnum]
                 , blocks[group==cc & acute.sc==ac & lab=="no AIDS mortality, gen het=0", jobnum])
    ## each transmission route
    inds[[2]] <- c(blocks[group==cc & acute.sc==ac & lab=="as fitted", jobnum]
                 , blocks[group==cc & acute.sc==ac & lab=="scale pre--couple, gen het=0" & bmb.sc %in% c(.1,10), jobnum])
    inds[[3]] <- c(blocks[group==cc & acute.sc==ac & lab=="as fitted", jobnum]
                 , blocks[group==cc & acute.sc==ac & lab=="scale extra--couple, gen het=0" & bme.sc %in% c(.1,10), jobnum])
    inds[[4]] <- c(blocks[group==cc & acute.sc==ac & lab=="as fitted", jobnum]
                 , blocks[group==cc & acute.sc==ac & lab=="scale within--couple, gen het=0" & bmp.sc %in% c(.1,10), jobnum])
    ## heterogeneity
    inds[[5]] <- c(blocks[group==cc & acute.sc==ac & lab=="as fitted", jobnum]
                 , blocks[group==cc & acute.sc==ac & lab=='gen cor=0' & het.gen.sd %in% 1:2, jobnum])
    ## assortativity (genetic, i.e. all route)
    inds[[6]] <- c(blocks[group==cc & acute.sc==ac & lab=='gen cor=0' & het.gen.sd==2, jobnum]
                 , blocks[group==cc & acute.sc==ac & lab=='gen cor=0.4' & het.gen.sd==2, jobnum]
                 , blocks[group==cc & acute.sc==ac & lab=='gen cor=0.8' & het.gen.sd==2, jobnum])
    inds
}
inds <- indsFxn(cc=15, ac = 1)
inds

mains <- c('A','B','C','D','E','F')
mains <- c('mortality','pre-couple \ncontact coefficient','extra-couple \ncontact coefficient',
           'intrinsic \ntransmission rate', 'heterogeneity \nin transmission', 'heterogeneity \nwith assortativity')

leg.cex <- .8
cc1 <- countries
cols <- gray(c(0,.4, .7))
pdf(file.path(dir.figs, 'Fig 2 - Counterfactual Summary.pdf'))
par(mfrow = c(2,3))
for(ii in 1:6) { ## each panel
    blocks[inds[[ii]], .(acute.sc, simj, jobnum, lab)] ## look at the meta data
    fls <- paste0(dir.results, "/Rdatas/CF-", formatC(inds[[ii]], width=6, flag=0), ".Rdata")
    ## load all files & combine
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
    tss[, cols:= cols[simNumplot]]
    tss[,ltys:= c(1, rep(2,10))[simNumplot]]
    tss[, plot(yr, sdp, type = 'n', lwd = 2, xlab = '', ylab = '', bty = 'n', xlim=c(1990, 2015), ylim = c(0,1), col = cols[jj],
                  main =mains[ii], las = 2)]
    tss[, lines(yr, sdp, col=cols[1], lty = ltys[1]), simNumplot]
    points(dframe.s$yr[dframe.s$group==cc1], dframe.s$psdc[dframe.s$group==cc1], pch = 19, col = 'black', cex = 1.1)     
    if(ii==1) legend('bottomleft', c('as fitted', 'no AIDS mortality'), lwd = 2:1, lty = 1, col = c('black', 'gray'), bty = 'n', cex = leg.cex)
    if(ii%in%2:3) legend('bottomleft', c('set to 0', 'as fitted', 'scaled X 10'), lwd = c(1,2,1), lty = c(1,1,2),
                         col = c( 'gray', 'black', 'gray'), bty = 'n', cex = leg.cex)
    if(ii%in%4) legend('bottomleft', c('scaled X 1/10', 'as fitted', 'scaled X 10'), lwd = c(1,2,1), lty = c(1,1,2),
                       col = c( 'gray', 'black',  'gray'), bty = 'n', cex = leg.cex)
    if(ii==5) legend('bottomleft', c('as fitted', 'std dev = 1', 'std dev = 2'), lwd = c(2,1,1), lty = c(1,1,2),
                     col = c('black',  'gray', 'gray'), bty = 'n', cex = leg.cex)
    if(ii==6) legend('bottomleft', c('ro = 0', 'ro = 0.4', 'ro = 0.8'), lwd = c(2,1,1), lty = c(1,1,2),
                     title = 'std dv = 2',
                     col = c('black',  'gray', 'gray'), bty = 'n', cex = leg.cex)
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
            sel <- which(cframe$group.ind==cc & cframe$simj %in% blocks$start[bltd[rr]]:blocks$end[bltd[rr]] & cframe$acute.sc==ac)
            xs <- cframe$bmb.sc[sel]
            xs[cframe$bmb.sc[sel]==0] <- .05 # since on a log scale
            lines(xs, (t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel],
                  col = cols[cc], type = 'l', pch = 19, lty = 1)
        }
        if(rr==3) {
            sel <- which(cframe$group.ind==cc & cframe$simj %in% blocks$start[bltd[rr]]:blocks$end[bltd[rr]] & cframe$acute.sc==ac)
            xs <- cframe$bme.sc[sel]
            xs[cframe$bme.sc[sel]==0] <- .05 # since on a log scale
            lines(xs, (t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel],
                  col = cols[cc], type = 'l', pch = 19, lty = 1)
        }
        if(rr %in% 4) {
            sel <- which(cframe$group.ind==cc & cframe$simj %in% (blocks$start[bltd[rr]]+1):blocks$end[bltd[rr]] & cframe$acute.sc==ac) 
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


####################################################################################################
## Figure SX for the manuscript - Heterogeneity and Assortativity
####################################################################################################
col.pl <- 'black'
ac <- 7 ## acute phase RH to use in Figure
for(one.file in c(F,T)) {
    leg.cex <- .57
    rmp <- colorRampPalette(c("yellow","red"))
    hazm <- c('bmb.sc','bme.sc','bmp.sc')   ## for each route get simjob with as fitted, 0, 10 scalars; acuteRH=7, no heterogeneity
    inds <- rbind(c(1, 91:92), c(1,91,119))
    bltd <- c(17, 24) ## blocks to show in summary figure (no AIDS mortality, 3 routes, within-couple, heterogeneity)
    mains <- paste0('(',LETTERS[1:2],')')
    mains <- c('transmission risk\nheterogeneity', 'transmission risk\nassortativity')
    tbl <- list(lab = blocks$lab[bltd], seq = list(inds[1,], inds[2,]))
    ltys <- 1:3
    lwds <- c(3,1.2,1.2)
    cx <- .9
    if(one.file) pdf(file.path(dir.figs,paste0('Figure 3 - HetAssort SummaryAc',ac,'.pdf')), w = 6.5, h = 3)
    for(cc1 in countries) {
        npan <- 2
        if(!one.file)  pdf(file.path(dir.figs,paste0('Figure X - HetAssort SummaryAc ', ds.nm[cc1],' Ac',ac,'.pdf')), w = 4, h = 3)
        layout(t(matrix(c(1:npan,npan*2+1,(npan+1):(2*npan),npan*2+1),npan+1,2)),w = c(rep(1,npan),.85))
        par(mar = c(3,1,2,0), oma = c(1,3,0,0), cex.lab = cx, cex.axis = cx, cex.main = cx, fg = col.pl, col.axis = col.pl,
            col.lab = col.pl, col = col.pl, col.main = col.pl)
        for(bb in 1:npan) {
            jst <- c(tbl$seq[[bb]])
            js <- which(cframe$simj %in% jst)
            js <- js[cframe$acute.sc[js]==ac & cframe$group.ind[js]==cc1] # select sims for this acute phase RH & country
            main <- mains[bb]
            js1 <- cframe$job[js[cframe$simj[js]==1]] # which line was as fitted? always simj=1
            set.labs(bb,js)
            yaxt <- ifelse(bb==1, T, F)
            cols <- c('black', 'gray', 'gray')
            ltys <- c(1,1,2)
            plot.sdp.nsub(js = js, leg = leg, js1 = js1, make.pdf = F, early.yr = 1985, show.pts = show.pts, pts.group = cc1,
                          main = main, cex.leg = .8, yaxt = yaxt, ylab = '', ltys = ltys, lwds = lwds, cols = cols,
                          title = legtitle, browse=F, col.pl = col.pl, show.leg = F, sep.leg = F)
            if(bb==1)   legend('bottomleft', c('as fitted', 'std dev = 1', 'std dev = 2'), lwd = c(2,1,1), lty = c(1,1,2),
                   col = c('black',  'gray', 'gray'), bty = 'n', cex = leg.cex)
            if(bb==2)   legend('bottomleft', c('st dev = 2', 'std dev = 2, corr = 0.4', 'std dev = 2, corr = 0.8'), lwd = c(2,1,1), lty = c(1,1,2),
                   col = c('black',  'gray', 'gray'), bty = 'n', cex = leg.cex)
        }
####################################################################################################
        ## Extract SDP's from 2008 for all scenarios to plot in 2nd row of figure
        par(mar = c(3,1,.5,0))
        cols <- rainbow(length(ds.nm))
        ylab <- '' #ifelse(rr==1,'SDP','')
        xmax <- c(1,10,10,10,2)
        for(rr in 1:npan) {
            yaxt <- ifelse(rr==1, T, F)
            logd <- ''
            if(rr==1)   xlim <- c(-.2,3.2) else xlim <- c(-.1, 1)
            plot(1,1, type = 'n', xlim = xlim, ylim = c(0,1), bty = 'n', axes = F,
                 xlab = '', ylab = ylab, main = '', log = logd)
            grcol <- rgb(t(col2rgb(gray(.6), alpha = .5)),max=255)
            segments(0,0,0,1, col = grcol, lwd = 3)
            if(yaxt) axis(2, at = seq(0,1,l=5), las = 2) else axis(2, at = seq(0,1,l=5), labels = NA)
            if(rr==1)   axis(1, 0:3) else axis(1, c(0,.4,.8))
            yrind <- which(t.arr[,2,1]==2008)
            for(cc in countries)   {
                if(rr == 1) {
                    sel <- which(cframe$group.ind==cc & cframe$simj %in% c(1,90:92) & cframe$acute.sc==ac)
                    lines(cframe$het.gen.sd[sel], (t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel],
                          col = cols[cc], type = 'l', pch = 19, lty = 1)
                }else{
                    sel <- which(cframe$group.ind==cc & cframe$simj %in% c(91,119,147) & cframe$acute.sc==ac)
                    lines(cframe$het.gen.cor[sel], (t.arr[yrind,'mm',sel] + t.arr[yrind,'ff',sel]) / t.arr[yrind,'inf.alive',sel],
                          col = cols[cc], type = 'l', pch = 19, lty = 1)
                }
            }
            if(rr==1)       mtext('standard deviation \nof risk distribution', side = 1, line = 3, adj = .5, cex = leg.cex, outer = F)
            if(rr==2)       mtext('inter-partner risk correlation', side = 1, line = 2, adj = .5, cex = leg.cex, outer = F)             
        }
        mtext('serodiscordant proportion (SDP)', side = 2, line = 2, adj = .5, cex = leg.cex, outer = T)
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

 
