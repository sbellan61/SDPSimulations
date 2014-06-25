####################################################################################################
## Makes a control file to send to the cluster with one line per model fitting assignment (each is a
## country-group / acute phase combination)
####################################################################################################

##  source('DHSFitMK.R')
setwd('/home1/02413/sbellan/DHSFitting/')    # set working directory
load('data files/ds.nm.all.Rdata')        # load country names
acute.sc <- c(1,5,7,10,25,30,40,50)     # acute phase relative hazards that will be used during DHS fits
#acute.sc <- c(7)
nac <- length(acute.sc)                 # how many?
countries <- 1:length(ds.nm)   # which countries to fit
countries <- countries[(ds.nm=='WA')]
#countries <- 14 #c(1,7,10:12)
ds.nm[countries]     # ones we're fitting
new.folder <- F     # always make a new folder for the result?
adapt <- F       # do an MCMC adaptive phase before main sampling?
sigma.ad <- T       # use sigma from a previously run adaptive phase (if not adapting here, otherwise use sigma from previous main phase)
basedir <- file.path('results','DHSFits')
if(!file.exists(basedir))        dir.create(basedir)

sink("WADHSFitControlFile.txt") ## write to Control File txt
for(cc in countries) { # for each country
    batchdirnm <- file.path(basedir, paste0(ds.nm[cc],'RealAcbackup')) ## create folder name
    if(!file.exists(batchdirnm))        dir.create(batchdirnm)
    if(!file.exists(paste(batchdirnm,'/routs/',sep='')))        dir.create(paste(batchdirnm,'/routs/',sep=''))
    for(ii in 1:nac) { ## for each acute phase RH
        jb <- ii                    # job num
        cmd <- paste("R CMD BATCH '--args jobnum=", jb,
                     ## country to analyze
                     " group.ind=", cc,
                     ## country to simulate
                     " s.epic=", cc, " s.demog=", cc, " s.bmb=", cc,
                     " s.bfb=", cc, " s.bme=", cc, " s.bfe=", cc, " s.bmp=", cc, " s.bfp=", cc,
                     ## Simulation parameters
                     " acute.sc=", acute.sc[ii], # acute phase to use in fitting
                     " bmp.sc=1 bfp.sc=1",       # this is fit to real data not simulation so these don't matter
                     " late.sc=1 aids.sc=1",  # no late or pre-aids death stages yet
                     " new.folder=", new.folder, " simul=FALSE maxN=1000 seed.bump=1", " sigma.ad=", sigma.ad,
                     " batchdirnm=\"", batchdirnm,"\"",
                     " psNonPar=TRUE infl.fac=50 death=TRUE bmb.sc=1 bfb.sc=1 bme.sc=1 bfe.sc=1 het.b=FALSE het.b.sd=0 het.b.cor=0 het.e=FALSE het.e.sd=0 het.e.cor=0 het.p=FALSE het.p.sd=0 het.p.cor=0 het.gen=FALSE het.gen.sd=0 het.gen.cor=0 het.beh=FALSE het.beh.sd=0 het.beh.cor=0 scale.by.sd=TRUE scale.adj=1 sample.tmar=FALSE tmar=(65*12):(111*12) each=200 tint=111*12 short.test=F all.cores=T nc=6",
                     ## fitting parameters
                     " adapt=", adapt,
                     " d.nburn=400 d.nthin=3 d.niter=1200 nburn=500 nthin=1 niter=7000 survive=T tell=100 low.coverage.arv=F partner.arv=F fsd.sens=F' DHSFitStarter.R ",
                     batchdirnm, '/routs/', ds.nm[cc],'RealAc',acute.sc[ii],'-j',jb,".Rout", sep = "")
        if(ii > 0) {
            cat(cmd)
            cat('\n')
          } } }
sink()
