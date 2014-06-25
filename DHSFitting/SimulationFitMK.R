####################################################################################################
## Makes a control file with each line sent to a cluster node to simulate DHS data & then fit it.
####################################################################################################
## source('SimulationFitMK.R')
setwd('/home1/02413/sbellan/DHSFitting/')    # set working directory
load('data files/ds.name.Rdata')        # load country names
## SIMULATIONS
countries <- 9                          # countries to do
nn <- 10                                # numer of simulations (all are going to be the same here)
acute.sc <- 7
acute.sc <- rep(acute.sc, nn)
group.ind <- countries
maxN <- rep(2000,nn)
seeds <-  (1:nn)*6
bmp.sc.vec <- 120 / (2*acute.sc + 118) # algebraic expression bp = (2/120)acute.sc * bp' + 118/120 * bp'
bfp.sc.vec <- 120 / (2*acute.sc + 118)

sink("SimulationFitControlFile.txt")
for(cc in countries) { # countries to do
    batchdirnm <- file.path('results','SimFits', paste0(ds.nm[cc],'SimAc',acute.sc[1],'N',maxN[1]))
    if(!file.exists(batchdirnm))        dir.create(batchdirnm)
    if(!file.exists(file.path(batchdirnm,'routs')))        dir.create(file.path(batchdirnm,'routs'))
    for(ii in 1:nn) { # for each run
        jb <- ii                    # job num
        cmd <- paste("R CMD BATCH '--args jobnum=", jb,
                     ## country to analyze
                     " group.ind=", cc,
                     ## country to simulate
                     " s.epic=", cc, " s.demog=", cc, " s.bmb=", cc,
                     " s.bfb=", cc, " s.bme=", cc, " s.bfe=", cc,
                     " s.bmp=", cc, " s.bfp=", cc,
                     ## Simulation parameters
                     " acute.sc=", acute.sc[ii], # acute phase
                     " bmp.sc=", bmp.sc.vec[ii], " bfp.sc=", bfp.sc.vec[ii],
                     " late.sc=1 aids.sc=1",  # no late or pre-aids death stages yet
                     " simul=TRUE",
                     " maxN=", maxN[ii],
                     " seed.bump=", seeds[ii],
                     " batchdirnm=\"", batchdirnm,"\"",
                     " psNonPar=TRUE infl.fac=50 death=TRUE bmb.sc=1 bfb.sc=1 bme.sc=1 bfe.sc=1 het.b=FALSE het.b.sd=0 het.b.cor=0 het.e=FALSE het.e.sd=0 het.e.cor=0 het.p=FALSE het.p.sd=0 het.p.cor=0 het.gen=FALSE het.gen.sd=0 het.gen.cor=0 het.beh=FALSE het.beh.sd=0 het.beh.cor=0 scale.by.sd=TRUE scale.adj=1 sample.tmar=FALSE tmar=(65*12):(111*12) each=200 tint=111*12 short.test=F all.cores=T nc=6",
                     ## fitting parameters
                     " adapt=T d.nburn=300 d.nthin=3 d.niter=700 nburn=500 nthin=1 niter=5000 survive=T tell=100 low.coverage.arv=F partner.arv=F fsd.sens=F' DHSFitStarter.R",
                     batchdirnm, '/routs/', ds.nm[cc],'Sim-',jb,".Rout", sep = "")
        if(ii > 0) {
            cat(cmd)
            cat('\n')
          } } }
sink()



