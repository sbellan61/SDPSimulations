rm(list=ls())
graphics.off()
setwd('/home1/02413/sbellan/DHSFitting')    # set working directory
lab <- 'EpiCurves'
wd <- file.path('results',lab)
outdir <- file.path('results',paste0(lab,'Summary'))
if(!file.exists(outdir)) dir.create(outdir)
load('data files/epic.Rdata')
load("data files/allDHSAIS.Rdata")         # DHS data
load("data files/ds.nm.all.Rdata") # country names

dat$wa <- dat$epic.nm %in%  c("Burkina Faso", "Cameroon", 'Cote dIvoire', "Ghana", "Guinea", "Liberia", "Mali", "Niger", "Senegal", "Sierra Leone")
countries <- unique(dat$epic.ind[!dat$wa])
countries.nm <- unique(dat$epic.nm[!dat$wa])
ncount <- length(countries)
rand <- order(rnorm(ncount))
cols <- rainbow(ncount)[rand]
for(fg in c('white','black')) {
  ## all prevalence ESA
  for(show.inf in c(T,F)) {
    pdf(file.path(outdir,paste0("ssa epic",' ART '[show.inf],fg,".pdf")), w = 4.5, h=3)
    ymax <- .3
    cex.nm <- .7
    par(mfrow = c(1,1), mar = c(.5,2,1.5,0), oma = c(2,2,0,0), fg = fg, col.lab = fg, col.axis = fg)
    plot(0, 0, xlim =c(1975,2030), pch = 19, bty = "n", ylim = c(0,ymax), type = "n",
         xlab = "", ylab = "prevalence", main = "", las = 2, xaxt="n", cex.axis = cex.nm)
    axis(1, seq(1975,2010, by = 5), NA)
    for(ii in 1:ncount) {
      ind <- countries[ii]
      show <- prev.all[,'cmc']/12 + 1900 > 1975 & prev.all[,'cmc']/12 + 1900 < 2011
                                        # lines(epicf[show,'cmc']/12+1900, epicf[show,ind], col = cols[ii], lty = 2)
      lines(prev.all[show,'cmc']/12+1900, prev.all[show,ind], col = cols[ii], lty = 1)
      if(show.inf)        lines(prev.inf[show,'cmc']/12+1900, prev.inf[show,ind], col = cols[ii], lty = 2)    
      text(2012, prev.all[112*12,ind], colnames(prev.all)[ind], pos = 4, cex = cex.nm, col = cols[ii])
    }
    axis(1, seq(1975,2010, by = 5), las = 2, cex.axis=cex.nm)
    mtext("HIV prevalence (15+)", side = 2, outer = T, line = 1, cex = cex.nm)
    dev.off()
  }
}

## all prevalence WA
countries <- unique(dat$epic.ind[dat$wa])
countries.nm <- unique(dat$epic.nm[dat$wa])
ncount <- length(countries)
rand <- order(rnorm(ncount))
cols <- rainbow(ncount)[rand]
for(fg in c('white','black')) {
  for(show.inf in c(T,F)) {
    pdf(file.path(outdir,paste0("wa epic",' ART '[show.inf],fg,".pdf")), w = 4.5, h=3)
    ymax <- .065
    cex.nm <- .7
    par(mfrow = c(1,1), mar = c(.5,2,1.5,0), oma = c(2,2,0,0), fg = fg, col.lab = fg, col.axis = fg)
    plot(0, 0, xlim =c(1975,2030), pch = 19, bty = "n", ylim = c(0,ymax), type = "n",
         xlab = "", ylab = "prevalence", main = "", las = 2, xaxt="n", cex.axis = cex.nm)
    axis(1, seq(1975,2010, by = 5), NA)
    for(ii in 1:ncount) {
      ind <- countries[ii]
      show <- prev.all[,'cmc']/12 + 1900 > 1975 & prev.all[,'cmc']/12 + 1900 < 2011
                                        # lines(epicf[show,'cmc']/12+1900, epicf[show,ind], col = cols[ii], lty = 2)
      lines(prev.all[show,'cmc']/12+1900, prev.all[show,ind], col = cols[ii], lty = 1)
      if(show.inf)        lines(prev.inf[show,'cmc']/12+1900, prev.inf[show,ind], col = cols[ii], lty = 2)    
      text(2012, prev.all[112*12,ind], colnames(prev.all)[ind], pos = 4, cex = cex.nm, col = cols[ii])
    }
    axis(1, seq(1975,2010, by = 5), las = 2, cex.axis=cex.nm)
    mtext("HIV prevalence (15+)", side = 2, outer = T, line = 1, cex = cex.nm)
    dev.off()
  }
}

## males
pdf(file.path(outdir,"ssa epic.pdf"), w = 4.5, h=4)
ymax <- .3
cex.nm <- .5
par(mfrow = c(2,1), mar = c(.5,2,1.5,0), oma = c(2,2,0,0))
plot(0, 0, xlim =c(1975,2011), pch = 19, bty = "n", ylim = c(0,ymax), type = "n",
     xlab = "", ylab = "prevalence", main = "", las = 2, xaxt="n", cex.axis = cex.nm)
axis(1, seq(1975,2010, by = 5), NA)
countries <- unique(dat$epic.ind[!dat$wa])
ncount <- length(countries)
cols <- rainbow(ncount)
for(ii in 1:ncount)
  {
    ind <- countries[ii]
    show <- epicm[,'cmc']/12 + 1900 > 1975 & epicm[,'cmc']/12 + 1900 < 2011
    lines(epicm[show,'cmc']/12+1900, epicm[show,ind], col = cols[ii], lty = 2)
    lines(epicm.all[show,'cmc']/12+1900, epicm.all[show,ind], col = cols[ii], lty = 1)
  }
mtext("males", side = 3, line = 0, cex = 1)
legend("topleft", colnames(epicm)[countries], col=cols, lwd=1, bty = "n", ncol = 2, cex = cex.nm)
plot(0, 0, xlim =c(1975,2011), pch = 19, bty = "n", ylim = c(0,ymax), type = "n",
     xlab = "", ylab = "prevalence", main = "", las = 2, xaxt="n", cex.axis = cex.nm)
axis(1, seq(1975,2010, by = 5), NA)
countries <- unique(dat$epic.ind[!dat$wa])
ncount <- length(countries)
cols <- rainbow(ncount)
for(ii in 1:ncount)
  {
    ind <- countries[ii]
    show <- epicf[,'cmc']/12 + 1900 > 1975 & epicf[,'cmc']/12 + 1900 < 2011
    lines(epicf[show,'cmc']/12+1900, epicf[show,ind], col = cols[ii], lty = 2)
    lines(epicf.all[show,'cmc']/12+1900, epicf.all[show,ind], col = cols[ii], lty = 1)
  }
mtext("females", side = 3, line = 0, cex = 1)
#legend("topleft", colnames(epicf)[countries], col=cols, lwd=1, bty = "n", ncol = 2, cex = cex.nm)
axis(1, seq(1975,2010, by = 5), las = 2, cex.axis=cex.nm)
mtext("population prevalence (15+)", side = 2, outer = T, line = 1, cex = cex.nm)
dev.off()


## make figure of epidemic curves
pdf(file.path(outdir,"wa epic.pdf"), w = 4.5, h=4)
ymax <- .08
cex.nm <- .5
par(mfrow = c(2,1), mar = c(.5,2,1.5,0), oma = c(2,2,0,0))
plot(0, 0, xlim =c(1975,2011), pch = 19, bty = "n", ylim = c(0,ymax), type = "n",
     xlab = "", ylab = "prevalence", main = "", las = 2, xaxt="n", cex.axis = cex.nm)
axis(1, seq(1975,2010, by = 5), NA)
countries <- unique(dat$epic.ind[dat$wa])
ncount <- length(countries)
cols <- rainbow(ncount)
for(ii in 1:ncount)
  {
    ind <- countries[ii]
    show <- epicm$cmc/12 + 1900 > 1975 & epicm$cmc/12 + 1900 < 2011
    lines(epicm$cmc[show]/12+1900, epicm[show,ind], col = cols[ii], lty = 2)
    lines(epicm.all$cmc[show]/12+1900, epicm.all[show,ind], col = cols[ii], lty = 1)
  }
mtext("males", side = 3, line = 0, cex = 1)
legend("topleft", names(epicm)[countries], col=cols, lwd=1, bty = "n", ncol = 2, cex = cex.nm)
plot(0, 0, xlim =c(1975,2011), pch = 19, bty = "n", ylim = c(0,ymax), type = "n",
     xlab = "year", ylab = "prevalence", main = "", las = 2, cex.axis = cex.nm)
countries <- unique(dat$epic.ind[dat$wa])
ncount <- length(countries)
cols <- rainbow(ncount)
for(ii in 1:ncount)
  {
    ind <- countries[ii]
    show <- epicf$cmc/12 + 1900 > 1975 & epicf$cmc/12 + 1900 < 2011
    lines(epicf$cmc[show]/12+1900, epicf[show,ind], col = cols[ii], lty = 2)
    lines(epicf.all$cmc[show]/12+1900, epicf.all[show,ind], col = cols[ii], lty = 1)
  }
mtext("females", side = 3, line = 0, cex = 1)
#legend("topleft", c("HIV prevalence","HIV prevalence * (1-ARV coverage)"), lty = 1:2, bty = "n", cex = cex.nm)
mtext("population prevalence (15+)", side = 2, outer = T, line = 1, cex = cex.nm)
dev.off()
