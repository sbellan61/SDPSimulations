####################################################################################################
## Collect, summarize, & visualize results from counterfactual simulations.
####################################################################################################
rm(list=ls())                           # clear workspace
gc()
library(plyr); library(mgcv)
## setwd('DHSProject/SDPSimulations')

source('PlotFunctions.R')                    # load functions to collect & plot results
source('SimulationFunctions.R')                   # load simulation functions
library(abind)                          # array binding
load('data files/ds.nm.all.Rdata')        # country names
load('data files/dframe.s.Rdata')        # SDP by survey
do.again <- F                           # collect results again (otherwise load cfs.Rdata)
show.pts <- T                           # show observed SDP in real DHS data

dir.results <- file.path('results','CounterAssort') # results locations
dir.results <- file.path('results','CounterAssortHiLo') # results locations
dir.figs <- file.path(dir.results, 'Figures')        # make a directory to store figures
if(!file.exists(dir.figs)) dir.create(dir.figs)      # create it
## load 'blocks' which gives info on all simulations within each country-acute group
load(file.path(dir.results, 'blocks.Rdata'))
## Load all results files in the results directory (all Rdata files except blocks & cfs)
fs <- list.files(pattern = '.Rdata', path = file.path(dir.results), recursive = T, full.names = T)
fs <- fs[!grepl('blocks',fs) & !grepl('cfs',fs) & !grepl('mardur',fs)]
#fs <- fs[grepl('Acute7',fs)]
print(length(fs))
## Order them by jobnumber
getjob <- function(x) as.numeric(sub('.Rdata', '', strsplit(x,c('-'))[[1]][3]))
ord <- order(sapply(fs, getjob))
fs <- fs[ord]
## coll: collects all results into data frame of input parameters, and array of time series
if(!file.exists(file.path(dir.results, 'cfs.Rdata')) | do.again) {
  cfs <- coll(fs, nc = 12, give.ev = F)
  save(cfs, file = file.path(dir.results, 'cfs.Rdata'))
}else{
  load(file.path(dir.results, 'cfs.Rdata'))
}
dim(cfs$cframe) ## dimensions & check for duplicate runs
print(paste(sum(duplicated(cfs$cframe$job)), 'duplicate jobs'))
cfs$cframe <- cfs$cframe[!duplicated(cfs$cframe$job),]
cfs$t.arr <-  cfs$t.arr[,,!duplicated(cfs$cframe$job)]
acutes <- unique(cfs$cframe$acute.sc)       # which acute phase relative hazards (RHs) were simulated
nac <- length(acutes)                   # how many
countries <- unique(cfs$cframe$group.ind)   # countries simulated
countries <- countries[order(countries)]
ngroup <- length(countries)             # how many
jtd <- which(!blocks$jobnum %in% cfs$cframe$job)    # jobs undone
print(paste("didn't do jobs:",paste(head(jtd,50), collapse=','))) # check to see if any jobs didn't complete
save(jtd, file='data files/AssortJobsToDo.Rdata')

####################################################################################################
## Couple Serostatus Summary Functions
sdpfx <- function(ser) sum(ser %in% 2:3) / sum(ser%in%1:3) ## serodiscordant proportion
sdp.lcifx <- function(ser)  ifelse(sum(ser%in%1:3)>0, unlist(binom.test(sum(ser%in%2:3),sum(ser%in%1:3))[4])[1], NA)
sdp.ucifx <- function(ser)  ifelse(sum(ser%in%1:3)>0, unlist(binom.test(sum(ser%in%2:3),sum(ser%in%1:3))[4])[2], NA)
prevfx <- function(ser) (sum(ser %in% 2:3) + 2*sum(ser==1)) / (2*length(ser)) ## prevalence
prev.lcifx <- function(ser)  ifelse(sum(ser%in%1:4), unlist(binom.test(sum(ser%in%1:3),sum(ser%in%1:4))[4])[1], NA)
prev.ucifx <- function(ser)  ifelse(sum(ser%in%1:4), unlist(binom.test(sum(ser%in%1:3),sum(ser%in%1:4))[4])[2], NA)
sdcfx <- function(ser) sum(ser %in% 2:3) / length(ser) ## prevalence of serodiscordance
cpcfx <- function(ser) sum(ser %in% 1) / length(ser)## prevalence of concordant positivity
 
## Extract these things from each simulation result as a function of time since couple formation
serostatus.by.mardur.mods <- function(fs, type = 'sdp') {
#  browser()
  load(fs)
  temp <- output$evout
  rm(output); gc()
  temp$is.sd <- temp$ser %in% 2:3
  temp$is.cp <- temp$ser %in% 1
  temp.pos <- temp[temp$ser %in% 1:3,]
  if(type=='sdp') mod <- gam(is.sd ~ s(mardur.mon, bs='cr', k = 5), family = binomial('logit'), temp.pos) ## cubic spline functions of SDP vs mardur
  if(type=='sdc') mod <- gam(is.sd ~ s(mardur.mon, bs='cr', k = 5), family = binomial('logit'), temp) ## serodiscordance prevalence vs mardur
  if(type=='cpc') mod <- gam(is.cp ~ s(mardur.mon, bs='cr', k = 5), family = binomial('logit'), temp) ## concordant positive prevalence vs mardur
  rm(temp, temp.pos); gc()
  return(list(mod))
}

wrp.serostatus.by.mardur.mods <- function(jobs, mc.cores = 12) { ## Parallelize it
  fss <- fs[jobs]
    out <- mclapply(fss, serostatus.by.mardur.mods, mc.cores = 12)
    names(out) <- jobs
    return(out)
  }

rmp <- colorRampPalette(c('black',"orange",'red',"purple"))
col.pl <- 'black'
ac <- 7 ## acute phase RH to use in Figure
blocks$phi <- blocks$phihi + blocks$phi.m

####################################################################################################
## Plot SDP in a survey year as a function of time since married
xs <- 0:400
for(sc in 1) { #unique(blocks$bmp.sc)) {
  pdf(file.path(dir.figs,paste0('Figure X - SDP vs mardur Ac',ac,' beta.sc', sc,'.pdf')), w = 6.5, h = 5)
  for(cc in 15) {#1:length(ds.nm)) { ## by country
    par(mfrow = c(2,2), mar = c(3,1,2,.5), oma = c(3,3,2,0), fg = col.pl, col.axis = col.pl, 'ps' = 12,
        col.lab = col.pl, col = col.pl, col.main = col.pl)
    for(hh in 1:4) { ## for each heterogeneity sd
      hsd <- unique(blocks$phi)[hh]
      js <- with(blocks, which(phi==hsd & group==cc & bmp.sc==sc)) ## simulations with that heterogeneity
      js1 <- with(blocks, which(phi==hsd & phihi==0 & group==cc & bmp.sc==sc)) ## the one with no correlation
      njs <- length(js)
      cols.rmp <- rmp(length(unique(blocks$rrhi.m)))
      cols <- cols.rmp[factor(blocks$rrhi.m[js])]
      ltys <- blocks$rel.assort[js]
      lwds <- rep(1,njs)
      main <- paste0('proportion high risk = ', hsd)
      js <- js[!js %in% jtd]
      temp.mods <- wrp.serostatus.by.mardur.mods(js)
      plot(0,0, type = 'n', xlab = '', ylab = '', main = main, xlim = c(0,25), ylim = c(0,1), las = 1, yaxt = 'n', bty = 'n')
      if(hh%%2==1) axis(2, at = seq(0,1, by = .2), las = 2)
      if(hh%%2==0) axis(2, at = seq(0,1, by = .2), lab = NA)
      if(hh==1) legend('bottomleft', leg = with(blocks[js,], paste0(rrhi.m,'X, hihi=', rel.assort,'-fold')), col = cols,
                                         lty = ltys, lwd = 1.5, cex = .8, title = '', ncol=4, bty = 'n')
      for(jj in 1:length(js)) { ## for each simulation, add lines to plot
        simj <- as.character(js[jj])
        ys <- predict(temp.mods[[simj]][[1]], data.frame(mardur.mon = xs, is.sd = NA), type = 'response')
        lines(xs/12, ys, col = cols[jj], lty = ltys[jj], lwd = 1.5)
      }
    }
    mtext('SDP', side = 2, line = 2, adj = .5, cex = 1, outer = T)
    mtext(paste0('couple duration (years) '), side = 1, line = 1, outer = T)
    mtext(ds.nm[cc], side = 3, line = 1, outer = T)
  }
  dev.off()
}

####################################################################################################
## Figure SDP time series by amount of intracouple correlation & intercouple heterogeneity
####################################################################################################
for(sc in unique(blocks$bmp.sc)) {
  pdf(file.path(dir.figs,paste0('Figure X - Assortativity Summary Ac',ac,' beta.sc', sc, '.pdf')), w = 4.5, h = 4)
  par(mfrow = c(2,2), mar = c(3,1,2,0), oma = c(1,3,2,0), fg = col.pl, col.axis = col.pl, 'ps' = 12,
      col.lab = col.pl, col = col.pl, col.main = col.pl)
  for(cc in 1:length(ds.nm)) { ## by country
    for(hh in 1:4) { ## for each amount of heterogeneity
      hsd <- c(0,1,2,3)[hh]
      js <- with(blocks, which(het.gen.sd==hsd & group==cc & bmp.sc==sc)) ## simulations with that heterogeneity
      js1 <- with(blocks, which(het.gen.sd==hsd & het.gen.cor==0 & group==cc & bmp.sc==sc)) ## the one with no correlation
      njs <- length(js)
      cols <- rmp(njs)
      ltys <- 1:4
      lwds <- rep(1,njs)
      main <- bquote(sigma[lambda]==.(hsd))
      yaxt <- hh%%2==1
      legtitle <- expression(rho[lambda])
      leg <- blocks$het.gen.cor[js]
      plot.sdp.nsub(js = js, leg = leg, js1 = js1, make.pdf = F, early.yr = 1985, show.pts = show.pts, pts.group = cc,
                    main = main, cex.leg = .8, yaxt = yaxt, ylab = '', ltys = ltys, lwds = lwds, cols = cols,
                    title = legtitle, col.pl = col.pl, show.leg = hh==1, sep.leg = F, browse=F)
    }
    mtext('SDP', side = 2, line = 2, adj = .5, cex = 1, outer = T)
    mtext(ds.nm[cc], side = 3, line = 1, outer = T)
  }
  graphics.off()
}
