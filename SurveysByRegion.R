######################################################################
#  Create a data frame with one row per country (dframe) or per
##  survey (dframe.s) with accompanying summary characteristics.
######################################################################
rm(list=ls())                           # clear workspace
library(plyr)
setwd('/home1/02413/sbellan/SDPSimulations/')     # setwd
#load("data files/ds.name.Rdata")        # country names
load("data files/ds.nm.all.Rdata") # country names
load("/home1/02413/sbellan/Rakai/SDPSimulations/data files/allDHSAIS.Rdata")         # DHS data
dato <- dat ## old data
load("data files/allDHSAIS.Rdata")         # DHS data
#load("data files/alldhs.Rdata")         # DHS data
load("data files/epic.Rdata")       # Infectious HIV prevalence in men
head(dat,2)
outdir <- file.path('results','PrevFigs')
if(!file.exists(outdir)) dir.create(outdir)
hazs <- c("bmb","bfb","bme","bfe","bmp","bfp")

####################################################################################################
## Get SDP by region
sdpfx <- function(ser) sum(ser %in% 2:3) / sum(ser%in%1:3)
sdp.lcifx <- function(ser)  ifelse(sum(ser%in%1:3)>0, unlist(binom.test(sum(ser%in%2:3),sum(ser%in%1:3))[4])[1], NA)
sdp.ucifx <- function(ser)  ifelse(sum(ser%in%1:3)>0, unlist(binom.test(sum(ser%in%2:3),sum(ser%in%1:3))[4])[2], NA)

prevfx <- function(ser) (sum(ser %in% 2:3) + 2*sum(ser==1)) / (2*length(ser))
prev.lcifx <- function(ser)  ifelse(sum(ser%in%1:4), unlist(binom.test(sum(ser%in%1:3),sum(ser%in%1:4))[4])[1], NA)
prev.ucifx <- function(ser)  ifelse(sum(ser%in%1:4), unlist(binom.test(sum(ser%in%1:3),sum(ser%in%1:4))[4])[2], NA)

rsdp <- ddply(dat, .(group,ds, region), summarise, sdp = sdpfx(ser), sdp.lci = sdp.lcifx(ser), sdp.uci = sdp.ucifx(ser),
              prev = prevfx(ser), prev.lci = prev.lcifx(ser), prev.uci = prev.ucifx(ser), yr = round(tint[1]/12+1900))
head(rsdp)

rsdp[rsdp$group=='Zambia',]

pdf(file.path(outdir, 'Regional SDP vs Prev.pdf'), w = 8, h = 6)
cols <- rainbow(nlevels(rsdp$group))
with(rsdp, {
  plot(prev, sdp, pch = 19, col = cols[as.numeric(group)],
                xlab = 'DHS prevalence', ylab = 'SDP')
  legend('topright', leg = levels(group), col = cols[1:nlevels(group)], pch = 19, ncol = 3)
} )
dev.off()

cis.s <- F
cis.p <- F
pdf(file.path(outdir, 'Regional by country SDP vs Prev.pdf'), w = 8, h = 6)
par(mfrow=c(4,4), mar = c(2,2,3,.5), oma = c(3,2,0,0))
for(cc in levels(rsdp$group)) {
  with(rsdp[rsdp$group==cc,], {
    dss <- unique(ds)
    plot(0,0, type = 'n', xlab = '', ylab = '', main = group[1], xlim = c(0,max(rsdp$prev)*1.1), ylim = c(0,1), las = 2)
    pchs <- 1:15
    regs <- factor(region)
    cols <- rainbow(length(dss), alpha = .45)
    cols[2] <- 'orange'
    cols.op <- rainbow(length(dss))
    cols.op[2] <- 'orange'
    for(ss in 1:length(dss)) {      ## Plot survey averages (what I used before)
      dst <- dss[ss]
      points(prevfx(dat$ser[dat$ds==dst]), sdpfx(dat$ser[dat$ds==dst]), pch = 19, col = cols.op[ss], cex = 2)
    }
    for(ss in 1:length(dss)) { ## Plot by region
      dst <- dss[ss]
      points(prev[ds==dst], sdp[ds==dst], col = cols[ss], pch = pchs[regs])
      if(cis.s) arrows(prev[ds==dst], sdp.lci[ds==dst], prev[ds==dst], sdp.uci[ds==dst], col = cols[ss], len = .02, angle = 90, code = 3)
      if(cis.p) arrows(prev.lci[ds==dst], sdp[ds==dst], prev.uci[ds==dst], sdp[ds==dst], col = cols[ss], len = .02, angle = 90, code = 3)      
      nans <- ds==dst & is.na(sdp) ## 0 prevalence regions
      if(sum(nans)>0) points(prev[nans], rep(0,sum(nans)), col = cols[ss], pch = pchs[regs[nans]])
    }
    legend('topright', leg=dss, pch = 19, col = cols, cex = .6)
  }
       )}
mtext('DHS prevalence', side = 1, line = 1, outer = T)
mtext('SDP', side = 2, line = 0, outer = T)
graphics.off()


####################################################################################################
## SDP by marital duration
pdf(file.path(outdir, 'SDP by cpl dur.pdf'), w = 8, h = 6)
par(mfrow=c(4,4), mar = c(2,2,3,.5), oma = c(3,2,0,0))
for(ii in 1:5) {
  ## yearly pooled marital duration
  breaks <- seq(0, 85, by = ii)
  labs <- breaks[-length(breaks)]
  dat$mdcn <- cut(dat$mardur.mon/12, breaks , lab = labs)
  dat$mdcn <- as.numeric(levels(dat$mdcn)[as.numeric(dat$mdcn)])
  rsdpm <- ddply(dat, .(group,ds, mdcn), summarise, sdp = sdpfx(ser), prev = prevfx(ser), yr = round(tint[1]/12+1900))
  rsdpm <- rsdpm[with(rsdpm, order(group, ds, mdcn)),]
  for(cc in levels(rsdp$group)) {
    with(rsdpm[rsdpm$group==cc,], {
      dss <- unique(ds)
      tds <- factor(ds)
      cols <- rainbow(length(dss))
      plot(0,0, type = 'n', xlab = '', ylab = '', main = group[1], xlim = c(0,25), ylim = c(0,1), las = 2)
      for(dst in dss)  {
#        if(cc=='Tanzania') browser()
        lines(mdcn[ds==dst], sdp[ds==dst], col = cols[as.numeric(tds[ds==dst])])
      }
      legend('topright', leg=dss, pch = 19, col = cols, cex = .6)
    }
         )}
  mtext(paste0('couple duration (',ii, ' year groupings)'), side = 1, line = 1, outer = T)
  mtext('SDP', side = 2, line = 0, outer = T)
}
graphics.off()

####################################################################################################
## SD & CC by marital duration
pdf(file.path(outdir, 'prop SDC & CCP by cpl dur.pdf'), w = 8, h = 6)
par(mfrow=c(4,4), mar = c(2,2,3,.5), oma = c(3,2,0,0))
for(ii in 1:5) {
  ## yearly pooled marital duration
  breaks <- seq(0, 85, by = ii)
  labs <- breaks[-length(breaks)]
  dat$mdcn <- cut(dat$mardur.mon/12, breaks , lab = labs)
  dat$mdcn <- as.numeric(levels(dat$mdcn)[as.numeric(dat$mdcn)])
  rsdpm <- ddply(dat, .(group,ds, mdcn), summarise, sdc = sum(ser %in% 2:3)/length(ser), ccp = sum(ser==1)/length(ser), yr = round(tint[1]/12+1900))
  rsdpm <- rsdpm[with(rsdpm, order(group, ds, mdcn)),]
  for(cc in levels(rsdp$group)) {
    with(rsdpm[rsdpm$group==cc,], {
      dss <- unique(ds)
      tds <- factor(ds)
      cols <- rainbow(length(dss))
      plot(0,0, type = 'n', xlab = '', ylab = '', main = group[1], xlim = c(0,25), ylim = c(0,.4), las = 2)
      for(dst in dss)  {
#        if(cc=='Tanzania') browser()
        lines(mdcn[ds==dst], sdc[ds==dst], col = cols[as.numeric(tds[ds==dst])])
        lines(mdcn[ds==dst], ccp[ds==dst], col = cols[as.numeric(tds[ds==dst])], lty = 2) 
      }
      legend('topright', leg=dss, pch = 19, col = cols, cex = .6)
      legend('topleft', leg=c("+-",'++'), lty = 1:2, cex = .6)
    }
         )}
  mtext(paste0('couple duration (',ii, ' year groupings)'), side = 1, line = 1, outer = T)
  mtext('proportion of couples by serostatus', side = 2, line = 1, outer = T)
}
graphics.off()


####################################################################################################
## look at same thing for some simulated data from each country Acute phase 7
for(ii in 1:nlevels(dat$group)) {
  cc <- levels(dat$group)[ii]
  print(paste0('adding simulated data from', cc))
  fls <- list.files(file.path('results','CounterFactual','Acute7',cc))
  fl <- fls[grepl(cc, fls)][1]
  load(file.path('results','CounterFactual','Acute7',cc,fl))
  if(ii==1)       ev <- output$ev else    ev <- rbind(ev, output$ev)
}
 
####################################################################################################
## SD & CC by marital duration BURUNDI SIM
pdf(file.path(outdir, 'SIMULATED prop SDC & CCP by cpl dur .pdf'), w = 8, h = 6)
par(mfrow=c(4,4), mar = c(2,2,3,.5), oma = c(3,2,2,0))
for(ii in 1:5) {
  ## yearly pooled marital duration
  breaks <- seq(0, 85, by = ii)
  labs <- breaks[-length(breaks)]
  ev$mdcn <- cut(ev$mardur.mon/12, breaks , lab = labs)
  ev$mdcn <- as.numeric(levels(ev$mdcn)[as.numeric(ev$mdcn)])
  rsdpm <- ddply(ev, .(group,ds, mdcn), summarise, sdc = sum(ser %in% 2:3)/length(ser), ccp = sum(ser==1)/length(ser), yr = round(tint[1]/12+1900))
  rsdpm <- rsdpm[with(rsdpm, order(group, ds, mdcn)),]
  for(cc in levels(rsdp$group)) {
    with(rsdpm[rsdpm$group==cc,], {
      dss <- unique(ds)
      tds <- factor(ds)
      cols <- rainbow(length(dss))
      plot(0,0, type = 'n', xlab = '', ylab = '', main = group[1], xlim = c(0,25), ylim = c(0,.6), las = 2)
      for(dst in dss)  {
#        if(cc=='Tanzania') browser()
        lines(mdcn[ds==dst], sdc[ds==dst], col = cols[as.numeric(tds[ds==dst])])
        lines(mdcn[ds==dst], ccp[ds==dst], col = cols[as.numeric(tds[ds==dst])], lty = 2) 
      }
      legend('topright', leg=dss, pch = 19, col = cols, cex = .6)
      legend('topleft', leg=c("+-",'++'), lty = 1:2, cex = .6)
    }
         )}
  mtext(paste0('couple duration (',ii, ' year groupings)'), side = 1, line = 1, outer = T)
  mtext('proportion of couples by serostatus', side = 2, line = 1, outer = T)
  mtext('SIMULATED DATA (as fit to DHS with acute 7)', side = 3, line = 0, outer = T)
}
graphics.off()


####################################################################################################
## SDP by marital duration SIMULATED
pdf(file.path(outdir, 'SIMULATED SDP by cpl dur.pdf'), w = 8, h = 6)
par(mfrow=c(4,4), mar = c(2,2,3,.5), oma = c(3,2,2,0))
for(ii in 1:5) {
  ## yearly pooled marital duration
  breaks <- seq(0, 85, by = ii)
  labs <- breaks[-length(breaks)]
  ev$mdcn <- cut(ev$mardur.mon/12, breaks , lab = labs)
  ev$mdcn <- as.numeric(levels(ev$mdcn)[as.numeric(ev$mdcn)])
  rsdpm <- ddply(ev, .(group,ds, mdcn), summarise, sdp = sdpfx(ser), prev = prevfx(ser), yr = round(tint[1]/12+1900))
  rsdpm <- rsdpm[with(rsdpm, order(group, ds, mdcn)),]
  for(cc in levels(rsdp$group)) {
    with(rsdpm[rsdpm$group==cc,], {
      dss <- unique(ds)
      tds <- factor(ds)
      cols <- rainbow(length(dss))
      plot(0,0, type = 'n', xlab = '', ylab = '', main = group[1], xlim = c(0,25), ylim = c(0,1), las = 2)
      for(dst in dss)  {
#        if(cc=='Tanzania') browser()
        lines(mdcn[ds==dst], sdp[ds==dst], col = cols[as.numeric(tds[ds==dst])])
      }
      legend('topright', leg=dss, pch = 19, col = cols, cex = .6)
    }
         )}
  mtext(paste0('couple duration (',ii, ' year groupings)'), side = 1, line = 1, outer = T)
  mtext('SDP', side = 2, line = 0, outer = T)
  mtext('SIMULATED DATA (as fit to DHS with acute 7)', side = 3, line = 0, outer = T)
}
graphics.off()
