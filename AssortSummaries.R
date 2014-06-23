####################################################################################################
## Collect, summarize, & visualize results from counterfactual simulations.
####################################################################################################
rm(list=ls())                           # clear workspace
library(plyr)

source('PlotFunctions.R')                    # load functions to collect & plot results
source('SimulationFunctions.R')                   # load simulation functions
library(abind)                          # array binding
load('data files/ds.nm.all.Rdata')        # country names
load('data files/dframe.s.Rdata')        # SDP by survey
do.again <- F                           # collect results again (otherwise load cfs.Rdata)
show.pts <- T                           # show observed SDP in real DHS data
## source('CounterFactualSummaries.R')

dir.results <- file.path('results','CounterAssort') # results locations
dir.figs <- file.path(dir.results, 'Figures')        # make a directory to store figures
if(!file.exists(dir.figs)) dir.create(dir.figs)      # create it
## load 'blocks' which gives info on all simulations within each country-acute group
load(file.path(dir.results, 'blocks.Rdata'))
## Load all results files in the results directory (all Rdata files except blocks & cfs)
fs <- list.files(pattern = '.Rdata', path = file.path(dir.results), recursive = T, full.names = T)
fs <- fs[!grepl('blocks',fs) & !grepl('cfs',fs)]
#fs <- fs[grepl('Acute7',fs)]
print(length(fs))
    
## coll: collects all results into data frame of input parameters, and array of time series
if(!file.exists(file.path(dir.results, 'cfs.Rdata')) | do.again) {
  cfs <- coll(fs, nc = 12, give.ev = T)
  cfs$cframe <- cfs$cframe[order(cfs$cframe$job),] # order by jobs
  cfs$t.arr <- cfs$t.arr[,,order(cfs$cframe$job)]          # ditto
  cfs$evout <- cfs$evout[,,order(cfs$cframe$job)]          # ditto  
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
nac <- length(acutes)                   # how many
countries <- unique(cframe$group.ind)   # countries simulated
countries <- countries[order(countries)]
ngroup <- length(countries)             # how many
jtd <- which(!blocks$jobnum %in% cframe$job)    # jobs undone
## print(paste("didn't do jobs:",paste(head(jtd,50), collapse=','))) # check to see if any jobs didn't complete
## save(jtd, file='data files/CFJobsToDo.Rdata')

cx <- .9
####################################################################################################
## Figure 2 for the manuscript
####################################################################################################
col.pl <- 'black'
ac <- 7 ## acute phase RH to use in Figure
## for(one.file in c(F,T)) {
rmp <- colorRampPalette(c("yellow","red"))
hazm <- c('bmb.sc','bme.sc','bmp.sc')   ## for each route get simjob with as fitted, 0, 10 scalars; acuteRH=7, no heterogeneity


pdf(file.path(dir.figs,paste0('Figure X - Assortativity Summary ', ds.nm[cc],' Ac',ac,'.pdf')), w = 4.5, h = 4)
par(mfrow = c(2,2), mar = c(3,1,2,0), oma = c(1,3,0,0), fg = col.pl, col.axis = col.pl, 'ps' = 12,
    col.lab = col.pl, col = col.pl, col.main = col.pl)
for(hh in 1:4) {
  hsd <- c(0,1,2,3)[hh]
  js <- with(blocks, which(het.gen.sd==hsd))# & het.gen.cor %in% c(0,.7,.9)))
  js1 <- which(blocks$het.gen.sd==hsd & blocks$het.gen.cor==0)
  njs <- length(js)
  cc <- blocks$group[1]
  cols <- rmp(njs)
  ltys <- 1:4
  lwds <- rep(1,njs)
  main <- bquote(sigma[lambda]==.(hsd))
  yaxt <- hh%%2==1
  legtitle <- 'intracouple correlation'
  leg <- blocks$het.gen.cor[js]
  plot.sdp.nsub(js = js, leg = leg, js1 = js1, make.pdf = F, early.yr = 1985, show.pts = show.pts, pts.group = cc,
                main = main, cex.leg = .8, yaxt = yaxt, ylab = '', ltys = ltys, lwds = lwds, cols = cols,
                title = legtitle, col.pl = col.pl, show.leg = hh==1, sep.leg = F, browse=F)
}
mtext('SDP', side = 2, line = 2, adj = .5, cex = 1, outer = T)
dev.off()


####################################################################################################
## Get SDP by region
sdpfx <- function(ser) sum(ser %in% 2:3) / sum(ser%in%1:3)
sdp.lcifx <- function(ser)  ifelse(sum(ser%in%1:3)>0, unlist(binom.test(sum(ser%in%2:3),sum(ser%in%1:3))[4])[1], NA)
sdp.ucifx <- function(ser)  ifelse(sum(ser%in%1:3)>0, unlist(binom.test(sum(ser%in%2:3),sum(ser%in%1:3))[4])[2], NA)

prevfx <- function(ser) (sum(ser %in% 2:3) + 2*sum(ser==1)) / (2*length(ser))
prev.lcifx <- function(ser)  ifelse(sum(ser%in%1:4), unlist(binom.test(sum(ser%in%1:3),sum(ser%in%1:4))[4])[1], NA)
prev.ucifx <- function(ser)  ifelse(sum(ser%in%1:4), unlist(binom.test(sum(ser%in%1:3),sum(ser%in%1:4))[4])[2], NA)

####################################################################################################
## Plot SDP in a survey year as a function of time since married
pdf(file.path(dir.figs,paste0('Figure X - SDP vs mardur ', ds.nm[cc],' Ac',ac,'.pdf')), w = 6.5, h = 4)
for(bb in 1:5) {
  breaks <- seq(0, 85, by = bb)
  labs <- breaks[-length(breaks)]
  par(mfrow = c(2,2), mar = c(3,1,2,0), oma = c(3,3,0,0), fg = col.pl, col.axis = col.pl, 'ps' = 12,
      col.lab = col.pl, col = col.pl, col.main = col.pl)
  for(hh in 1:4) {
    hsd <- c(0,1,2,3)[hh]
    js <- with(blocks, which(het.gen.sd==hsd))# & het.gen.cor %in% c(0,.7,.9)))
    js1 <- which(blocks$het.gen.sd==hsd & blocks$het.gen.cor==0)
    njs <- length(js)
    cc <- blocks$group[1]
    cols <- rmp(njs)
    cols[1] <- 'black'
    ltys <- 1:4
    lwds <- rep(1,njs)
    main <- bquote(sigma[lambda]==.(hsd))
    yaxt <- ifelse(hh%%2==1, 's','n')
    plot(0,0, type = 'n', xlab = '', ylab = '', main = main, xlim = c(0,25), ylim = c(0,1), las = 2, yaxt = yaxt, bty = 'n')
    if(hh==1) legend('topright', leg = blocks$het.gen.cor[js], col = cols, lty = ltys, cex = .6, title = expression(ro), ncol=4)
    for(jj in 1:length(js)) {
      temp <- cfs$e.arr[,,js[jj]]
      te <- cut(temp[,'mardur.mon']/12, breaks , lab = labs)
      temp <- cbind(temp, mdcn = as.numeric(levels(te)[te]))
      rsdpm <- ddply(as.data.frame(temp), .(mdcn), summarise, sdp = sdpfx(ser), prev = prevfx(ser), yr = round(tint[1]/12+1900))
      with(rsdpm, lines(mdcn, sdp, col = cols[jj], lty = ltys[jj]))
    }
  }
    mtext('SDP', side = 2, line = 2, adj = .5, cex = 1, outer = T)
    mtext(paste0('couple duration (',bb, ' year groupings)'), side = 1, line = 1, outer = T)
}
graphics.off()
