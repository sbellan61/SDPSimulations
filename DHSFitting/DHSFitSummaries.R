#################################################################################################### 
##  Collect results from fitting all DHS data sets across assume that
##  acute relative hazards of 1-50.
####################################################################################################
## Steve Bellan, 2013
####################################################################################################
rm(list=ls())
setwd('/home1/02413/sbellan/DHSFitting')    # set working directory
lab <- 'DHSFits'
wd <- file.path('results',lab)
outdir <- file.path('results',paste0(lab,'Summary'))
if(!file.exists(outdir)) dir.create(outdir)
fls <- list.files(wd, full.names=T)
## load a preliminary workspace to get dimensions of output for creating new arrays
dirs <- list.files(fls[1], full.names=T)
dirs <- dirs[!grepl('routs', dirs)]     # ignore Rout output files
load(file.path(dirs[1],'workspace.Rdata'))
library(coda)

## use old WA directory for now
fls <- fls[!grepl('shortchain', fls)]
#fls[grepl('WA',fls)] <- 'results/DHSFitsold/WARealAc/'
## Collect all estimated of proportion of transmission by routes, and
## route-gender specific prevalence standardized hazards into arrays
## for sensitivity analysis by acute phase by country
fls                                  
for(fff in 1:length(fls)) { # looping through these countries
    print(paste('extracting results from', fls[fff]))
        dirs <- list.files(fls[fff],full.name=T) # names of folders in each country folder (each is a different acute phase relative hazard)
        dirs <- dirs[!grepl('routs', dirs)]
#        if(fff==14) dirs <- dirs[2] # for WA which hasn't finished yet use old sim
        for(ddd in 1:length(dirs)) {
              if(fff == 1 & ddd==1)     # if first iteration of loop, set up array structure
                { ## out.arr.t stores median and credible intervals of traced parameters
                  out.arr.t <- array(NA, dim = c(nrow(pars)+12, ncol(pars), 8, length(fls))) # + 8 is for added gender geometric means & contact mixing pars below
                  colnames(out.arr.t) <- colnames(pars)
                  in.arr.t <- array(NA, dim = c(8, length(fls),3)) # in.arr.t stores country, acute relative hazard, and muttivariate gelman-rubin diagnostic
                }
              if(sum(grepl('Figure 6', list.files(dirs[ddd])))>0) {
                      load(file.path(dirs[ddd],'workspace.Rdata'))
                      in.arr.t[ddd,fff,1] <- group # country
                      in.arr.t[ddd,fff,2] <- acute.sc # acute relative hazard
                      in.arr.t[ddd,fff,3] <- gelout[[2]]
                      ## get geometric mean parameters for each transmission route and contact rates
                      ## (pre-/extra- transmission coefficients / within-couple hazard)
                      bmb.vec <- unlist(mcmc.out[,"bmb"]) # extract MCMC chains from list into one long vector
                      bfb.vec <- unlist(mcmc.out[,"bfb"])
                      bmp.vec <- unlist(mcmc.out[,"bmp"])
                      bfp.vec <- unlist(mcmc.out[,"bfp"])    
                      bme.vec <- unlist(mcmc.out[,"bme"])
                      bfe.vec <- unlist(mcmc.out[,"bfe"])
                      bb <- sqrt(bmb.vec*bfb.vec) # pre- geometric mean
                      be <- sqrt(bme.vec*bfe.vec) # extra- geometric mean
                      bp <- sqrt(bmp.vec*bfp.vec) # within- geometric mean
                      rr.ep <- be / bp            # extra- contact mixing (of geometric means)
                      rr.bp <- bb / bp            # pre- contact mixing (of geometric means)
                      rr.bp.m <- bmb.vec / bp     # pre- contact mixing (male) (using geometric mean for within-rate)
                      rr.bp.f <- bfb.vec / bp     # pre- contact mixing (female) (using geometric mean for within-rate)
                      rr.eb <- be / bb            # extra- / pre-
                      cmb <- bmb.vec / bmp.vec    # male pre- contact
                      cfb <- bfb.vec / bfp.vec    # female pre- contact
                      cme <- bme.vec / bmp.vec    # male extra- contact
                      cfe <- bfe.vec / bfp.vec    # female extra- contact
                      gmpars <- data.frame(cmb=cmb, cfb=cfb, cme=cme, cfe= cfe, bb=bb,be=be,bp=bp, rr.ep = rr.ep, rr.bp = rr.bp, rr.eb = rr.eb, rr.bp.m, rr.bp.f)
                      gmpars <- t(apply(gmpars, 2, function(x) quantile(x, c(.025,.5,.975)))) # get median/credible intervals
                      pars <- rbind(gmpars, pars) # add to pars
                      out.arr.t[,,ddd,fff] <- pars    # estimated hazards & route-specific contributions to infections
                      if(fff==1 & ddd==1)       rownames(out.arr.t) <- rownames(pars) # set rownames
                  }
            }
    }

in.arr <- in.arr.t ## renamed above because load(workspace) loaded old version of in.arr
out.arr <- out.arr.t

## reorder everything by acute phase infectivity (it was ordered by character order not numeric)
ord <- order(as.numeric(in.arr[,1,2]))
in.arr <- in.arr[ord,,]
out.arr <- out.arr[,,ord,]


## what does it look like?
#in.arr

in.arr[,,1]
in.arr[,,2]
in.arr[,,3]
apply(in.arr[,,3],2,max)
out.arr[hazs,2,,14]
ds.nm[is.na(in.arr[1,,1])]           # countries left
## save this for simulation later
save(out.arr, in.arr, file = "data files/pars.arr.ac.Rdata")
save(out.arr, in.arr, file = "/home1/02413/sbellan/SDPSimulations/data files/pars.arr.ac.Rdata")

load(file = "/home1/02413/sbellan/SDPSimulations/data files/pars.arr.ac.Rdata")

hazs <- c("bmb","bfb","bme","bfe","bmp","bfp")
out.arr[hazs,,,] <- out.arr[hazs,,,] * 12 * 100 # per 100 person years
## what's the distribution of Gelman-Rubin Diagnostics?
pdf(file.path(outdir,'gel-rubin hist.pdf'))
hist(as.numeric(in.arr[,,3]), xlab = "Gelman-Rubin Multivariate Diagnostic", col = 'black', main='Fitting all countries across acute relative hazard 1-50')
dev.off()

## ## for plotting, make hazards in units of per 100 person-years
## hazs <- c("bmb","bfb","bme","bfe","bmp","bfp")
## xmax <- 50 #max(as.numeric(in.arr[,2,2]), na.rm=T)
## cols <- rainbow(ncol(in.arr))
## ##  Plot fitted to hazards for all gender route specific transmission routes across the range of
## ##  assumed acute relative hazards.
## for(arws in c(T,F)) {
##   pdf(file.path(outdir,paste0('hazs vs acute','arrows'[arws],'.pdf')), w = 4.5, h = 6)
##                                         #par(oma = c(0,0,3,0), mar = c(4,4,2,0))
##   ylabs <- rep(paste(c('pre-couple', 'extra-couple', 'within-couple'),'beta'), each = 2)
##   for(zz in c(T,F)) { ## zz=T  make all y-axis limits the same, =F zoom in
##     par(oma = c(0,0,3,0), mar = c(4,4,2,0))
##     layout(matrix(c(1,3,5,2,4,6,7,7,7), nr = 3, nc = 3), widths = c(1,1,.5), heights = rep(1,3))
##     for(ii in 1:length(hazs)) {     # for each route
##       hh <- hazs[ii]
##       if(zz) ymax.ind <- hazs else ymax.ind <-  hh
##       if(zz) ymax.ind <- hazs else ymax.ind <-  hh
##       if(arws) ymax <- max(out.arr[ymax.ind,,,],na.rm=T) else ymax <- max(out.arr[ymax.ind,2,,],na.rm=T)
##       plot(0,0, type = 'n', bty = 'n', xlab = 'acute RR', ylab=ylabs[ii], ylim = c(0, ymax),
##            xlim = c(0, xmax))
##       for(fff in 1:dim(out.arr)[4]) { # for each country
##         sel <- which(!is.na(out.arr[hh,2,,fff]))
##         xs <- jitter(as.numeric(in.arr[sel,fff,2]), a = 1)
##         lines(xs, out.arr[hh,2,sel,fff],, col = cols[fff]) # connect medians with a line
##         for(ddd in sel) { # add error bars for  credible intervals at each median estimate
##           if(arws) {
##             arrows(xs, out.arr[hh,1,ddd,fff],
##                    xs, out.arr[hh,3,ddd,fff],
##                    code = 3, len = .02, angle = 90, col = cols[fff])
##           }
##         }
##       }
##     }
##     mtext('women', side = 3, line = -.5, outer = T, adj = .8)
##     mtext('men', side = 3, line = -.5, outer = T, adj = .25)
##     mtext('prevalence-standardized hazards (per 100 person-years)', side = 3, line = 1.5, outer = T)
##     par(mar=c(0,1,0,0))
##     plot(0,0, type = 'n', axes = F, bty = 'n')
##     legend('topleft', in.arr[1,,1], pch = 19, col = cols, cex = .7, ncol = 1, bty = 'n')
##   }
##   dev.off()       
## }


## for plotting, show contact coefficients & intrinsic transmisison rates in units of per 100 person-years
hazs <- c("bmb","bfb","bme","bfe","bmp","bfp")
xmax <- 50 #max(as.numeric(in.arr[,2,2]), na.rm=T)
cols <- rainbow(ncol(in.arr))
##  Plot fitted to hazards for all gender route specific transmission routes across the range of
##  assumed acute relative hazards.
for(arws in c(T,F)) {
  pdf(file.path(outdir,paste0('betas vs acute','arrows'[arws],'.pdf')), w = 4.5, h = 6)
                                        #par(oma = c(0,0,3,0), mar = c(4,4,2,0))
  ylabs <- c(expression(beta['M,b']),expression(beta['F,b']), expression(beta['M,e']), expression(beta['F,e']), expression(beta['M,p']),   expression(beta['F,p']))
  for(zz in c(F)) { ## zz=T  make all y-axis limits the same, =F zoom in
    par(oma = c(2,.3,1,0), mar = c(2,4,.5,1), xpd=NA)
    layout(matrix(c(1,3,5,2,4,6,7,7,7), nr = 3, nc = 3), widths = c(1,1,.5), heights = rep(1,3))
    for(ii in 1:length(hazs)) {     # for each route
      hh <- hazs[ii]
      if(zz) ymax.ind <- hazs else ymax.ind <-  hh
      if(arws) ymax <- max(out.arr[ymax.ind,,,],na.rm=T) else ymax <- max(out.arr[ymax.ind,2,,],na.rm=T)
      par(xpd=NA)
      ylab <- ylabs[ii]
      plot(0,0, type = 'n', bty = 'n', xlab = '', ylab=ylab, ylim = c(0, ymax),
           xlim = c(0, xmax), las = 2, axes = F)
      par(xpd=T)
      if(ii %in% 5:6) axis(1, at = seq(0,50, by = 10)) else  axis(1, at = seq(0,50, by = 10), label = NA)
      axis(2, at = pretty(c(0,ymax*.8),5)) 
      for(fff in 1:dim(out.arr)[4]) { # for each country
        sel <- which(!is.na(out.arr[hh,2,,fff]))
        xs <- jitter(as.numeric(in.arr[sel,fff,2]), a = 1)
        lines(xs, out.arr[hh,2,sel,fff],, col = cols[fff]) # connect medians with a line
        for(ddd in sel) { # add error bars for  credible intervals at each median estimate
          if(arws) {
            arrows(xs, out.arr[hh,1,ddd,fff],
                   xs, out.arr[hh,3,ddd,fff],
                   code = 3, len = .02, angle = 90, col = cols[fff])
          }
        }
      }
    }
    mtext('female', side = 3, line = -.5, outer = T, adj = .7)
    mtext('male', side = 3, line = -.5, outer = T, adj = .23)
    mtext('assumed acute phase to chronic phase relative hazard', side = 1, line = .8, outer = T, adj = .5, cex = .7)
    #mtext('prevalence-standardized hazards (per 100 person-years)', side = 3, line = 1.5, outer = T)
    par(mar=c(0,1,0,0))
    plot(0,0, type = 'n', axes = F, bty = 'n')
    legend('topleft', in.arr[1,,1], pch = 19, col = cols, cex = .8, ncol = 1, bty = 'n')
  }
  dev.off()       
}


## for plotting, show contact coefficients & intrinsic transmisison rates in units of per 100 person-years
hazs <- c("cmb","cfb","cme","cfe","bmp","bfp")
xmax <- 50 #max(as.numeric(in.arr[,2,2]), na.rm=T)
cols <- rainbow(ncol(in.arr))
##  Plot fitted to hazards for all gender route specific transmission routes across the range of
##  assumed acute relative hazards.
for(arws in c(T,F)) {
  pdf(file.path(outdir,paste0('cs and btas vs acute','arrows'[arws],'.pdf')), w = 4.5, h = 6)
                                        #par(oma = c(0,0,3,0), mar = c(4,4,2,0))
  ylabs <- c(expression(c['M,b']),expression(c['F,b']), expression(c['M,e']), expression(c['F,e']), expression(beta['M,p']),   expression(beta['F,p']))
  #ylabs <- rep(c(expression(c['pre-couple']), expression(c['extra-couple']), expression(beta^'*')),each=2)
  for(zz in c(F)) { ## zz=T  make all y-axis limits the same, =F zoom in
    par(oma = c(2,.3,1,0), mar = c(2,4,.5,1), xpd=NA)
    layout(matrix(c(1,3,5,2,4,6,7,7,7), nr = 3, nc = 3), widths = c(1,1,.5), heights = rep(1,3))
    for(ii in 1:length(hazs)) {     # for each route
      hh <- hazs[ii]
      if(zz) ymax.ind <- hazs else ymax.ind <-  hh
      if(arws) ymax <- max(out.arr[ymax.ind,,,],na.rm=T) else ymax <- max(out.arr[ymax.ind,2,,],na.rm=T)
      #ylab <- ifelse(ii %in% c(1,3,5), ylabs[ii],'')
      ylab <- ylabs[ii]
      par(xpd=NA)
      plot(0,0, type = 'n', bty = 'n', xlab = '', ylab=ylab, ylim = c(0, ymax),
           xlim = c(0, xmax), las = 2, axes = F)
      if(ii %in% 5:6) axis(1, at = seq(0,50, by = 10)) else  axis(1, at = seq(0,50, by = 10), label = NA)
      axis(2, at = pretty(c(0,ymax*.8),5))
      par(xpd=T)
      for(fff in 1:dim(out.arr)[4]) { # for each country
        sel <- which(!is.na(out.arr[hh,2,,fff]))
        xs <- jitter(as.numeric(in.arr[sel,fff,2]), a = 1)
        lines(xs, out.arr[hh,2,sel,fff],, col = cols[fff]) # connect medians with a line
        for(ddd in sel) { # add error bars for  credible intervals at each median estimate
          if(arws) {
            arrows(xs, out.arr[hh,1,ddd,fff],
                   xs, out.arr[hh,3,ddd,fff],
                   code = 3, len = .02, angle = 90, col = cols[fff])
          }
        }
      }
    }
    mtext('female', side = 3, line = -.5, outer = T, adj = .7)
    mtext('male', side = 3, line = -.5, outer = T, adj = .23)
    mtext('assumed acute phase to chronic phase relative hazard', side = 1, line = .8, outer = T, adj = .5, cex = .7)
    #mtext('prevalence-standardized hazards (per 100 person-years)', side = 3, line = 1.5, outer = T)
    par(mar=c(0,1,0,0))
    plot(0,0, type = 'n', axes = F, bty = 'n')
    legend('topleft', in.arr[1,,1], pch = 19, col = cols, cex = .8, ncol = 1, bty = 'n')
  }
  dev.off()       
}

##  Plot fitted proportion of transmission having occurred from each gender-specific route across
##  the range of assumed acute relative hazards.
piUs <- grepl('pi', rownames(pars)) & grepl('U', rownames(pars)) # get names of fitted proportions
piUs <- rownames(pars)[piUs]
ylabs <- rep(c('pre-couple', 'extra-couple', 'within-couple', 'within-couple ACUTE', 'within-couple CHRONIC'),2)
for(arws in c(F,T)) {
  pdf(file.path(outdir,paste0('pis vs acute','arrows'[arws],'.pdf')), w = 8, h = 5)
  par(mfrow = c(2,5), mar = c(3,4.5,.5,.5), oma = c(0,2,3,0))
  for(ii in 1:length(piUs)) {             # for each route
    hh <- piUs[ii]
    plot(0,0, type = 'n', bty = 'n', xlab = '', ylab=ylabs[ii], ylim = c(0,1), xlim = c(0, xmax), xaxt = 'n')
    axis(1, seq(0,50,b=5), label = NA)
    axis(1, c(0,25,50))
    for(fff in 1:dim(out.arr)[4]) { # for each country plot a line connecting the medians at each acute relative hazard
      lines(as.numeric(in.arr[,fff,2]), out.arr[hh,2,,fff], col = cols[fff])
      for(ddd in 1:dim(out.arr)[3]) { # and also plot credible interval error bars
        xs <- as.numeric(in.arr[ddd,fff,2])
        if(arws) {
          arrows(xs, out.arr[hh,1,ddd,fff],
                 xs, out.arr[hh,3,ddd,fff],
                 code = 3, len = .02, angle = 90, col = cols[fff])
        }
      }
    }
  }
  legend('top', in.arr[1,,1], pch = 19, col = cols, cex = .5, ncol = 2)
  mtext('men', side = 2, line = 0.2, outer = T, adj = .79)
  mtext('women', side = 2, line = 0.2, outer = T, adj = .25)
  mtext('proportion of transmission to each sex by route', side = 3, line = 1, outer = T)
  mtext('assumed acute phase to chronic phase relative hazard', side = 1, line = -1.2, outer = T, adj = .5, cex = .7)
  dev.off()       
}

####################################################################################################
## Create estimates of transmission rates by country & acute phase
hpars <- out.arr[c('bmp','bfp','bp'),,,]
hpars['bp',,,] <- 12*100*hpars['bp',,,] ## per 100 person years (other two rows have already been scaled)
hpars <- signif(hpars,2)
#hpars
paste0(hpars[,2,,], ' (', hpars[,1,,], ', ', hpars[,3,,], ')')
for(cc in 1:length(ds.nm)) { # countries
  for(pp in 1:3) { # bmp, bfp, bp
      if(cc==1 & pp == 1) {
        tab <- c(ds.nm[cc], c('male','female','mean')[pp], paste0(hpars[pp,2,,cc], ' (', hpars[pp,1,,cc], ', ', hpars[pp,3,,cc], ')'))
      }else{
        tab <- rbind(tab, c(ds.nm[cc], c('male','female','mean')[pp], paste0(hpars[pp,2,,cc], ' (', hpars[pp,1,,cc], ', ', hpars[pp,3,,cc], ')')))
      }
    }
}
tab[1:nrow(tab) %% 3 !=1,1] <- ''
colnames(tab) <- c('country', 'parameter', in.arr[,1,2])
tab
write.csv(tab, file=file.path(outdir,'transmission rate table.csv'))


####################################################################################################
## Create estimates of transmission rates by country for acute phase = 7
hpars <- 12*100*out.arr[c('bmp','bfp','bp'),,,]
#hpars['bp',,,] <- 12*100*hpars['bp',,,] ## per 100 person years (other two rows have already been scaled)
hpars <- signif(hpars,2)
hpars <- hpars[,,3,]
#hpars
for(cc in 1:length(ds.nm)) { # countries
  for(pp in 1:3) { # bmp, bfp, bp
      if(cc==1 & pp == 1) {
        tab <- c(ds.nm[cc], c('male','female','mean')[pp], paste0(hpars[pp,2,cc], ' (', hpars[pp,1,cc], ', ', hpars[pp,3,cc], ')'))
      }else{
        tab <- rbind(tab, c(ds.nm[cc], c('male','female','mean')[pp], paste0(hpars[pp,2,cc], ' (', hpars[pp,1,cc], ', ', hpars[pp,3,cc], ')')))
      }
    }
}
tab[1:nrow(tab) %% 3 !=1,1] <- ''
colnames(tab) <- c('country', 'parameter', in.arr[,1,2])
tab

write.csv(tab, file=file.path(outdir,'transmission rate table Ac7.csv'))

####################################################################################################
## barplot

##  Plot fitted proportion of transmission having occurred from each gender-specific route across
##  the range of assumed acute relative hazards.
ac <- 7
ac.ind <- which(in.arr[,1,2]==ac)
beside <- F
show <- c("pibUA","pieUA","pipUA")
mbdn.all <- out.arr[show,,ac.ind,]
mbdn.all[,2,] <- apply(mbdn.all[,2,], 2, function(x) {x/sum(x)}) # normalize medians to add to 1 (don't have to b/c they come from accross chains)
show <- c("piUbA","piUeA","piUpA")
fbdn.all <- out.arr[show,,ac.ind,]
fbdn.all[,2,] <- apply(fbdn.all[,2,], 2, function(x) {x/sum(x)})
space <- .2
if(beside) space <- c(0,.5)
ns <- length(ds.nm)
pdf(file.path(outdir,paste0('pisAc',ac,'.pdf')), w = 8, h = 5)
par(mar=c(4,6,3,2), mfrow = c(1,2), oma = c(0,0,1,0))
cols <- c("dark gray","red", "blue")
xlim <- c(0,1)
bp.temp <- barplot(mbdn.all[,2,ns:1], beside = beside, names.arg = ds.nm[ns:1], las =2, space = space,
                   horiz = T, col = cols, xlim = xlim, main = "males")
bp.temp <- barplot(fbdn.all[,2,ns:1], beside = beside, names.arg = ds.nm[ns:1], las =2, space = space,
                   horiz = T, col = cols, xlim = xlim, main = "females")
mtext("proportion of infected individuals", side = 1, line = -1.3, outer=T,
      cex = 1.3)
dev.off()


  par(mfrow = c(2,5), mar = c(3,4.5,.5,.5), oma = c(0,2,3,0))
  for(ii in 1:length(piUs)) {             # for each route
    hh <- piUs[ii]
    plot(0,0, type = 'n', bty = 'n', xlab = '', ylab=ylabs[ii], ylim = c(0,1), xlim = c(0, xmax), xaxt = 'n')
    axis(1, seq(0,50,b=5), label = NA)
    axis(1, c(0,25,50))
    for(fff in 1:dim(out.arr)[4]) { # for each country plot a line connecting the medians at each acute relative hazard
      lines(as.numeric(in.arr[,fff,2]), out.arr[hh,2,,fff], col = cols[fff])
      for(ddd in 1:dim(out.arr)[3]) { # and also plot credible interval error bars
        xs <- as.numeric(in.arr[ddd,fff,2])
        if(arws) {
          arrows(xs, out.arr[hh,1,ddd,fff],
                 xs, out.arr[hh,3,ddd,fff],
                 code = 3, len = .02, angle = 90, col = cols[fff])
        }
      }
    }
  }
  legend('top', in.arr[1,,1], pch = 19, col = cols, cex = .5, ncol = 2)
  mtext('men', side = 2, line = 0.2, outer = T, adj = .79)
  mtext('women', side = 2, line = 0.2, outer = T, adj = .25)
  mtext('proportion of transmission to each sex by route', side = 3, line = 1, outer = T)
  mtext('assumed acute phase to chronic phase relative hazard', side = 1, line = -1.2, outer = T, adj = .5, cex = .7)
  dev.off()       
}
