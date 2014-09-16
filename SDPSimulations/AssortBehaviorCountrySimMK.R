####################################################################################################
## Makes control files for each analysis within which each line giving one R CMD BATCH command line
## to run on a cluster.
####################################################################################################
rm(list=ls())                                  # clear workspace
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/DHSProject/SDPSimulations/')
load("data files/ds.nm.all.Rdata") # country names
load('data files/pars.arr.ac.Rdata')    # load acute phase relative hazards used to fit (in.arr[,,2])
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') #  transmission coefficient names, for convenience
nc <- 12                                       # core per simulation
## source('AssortBehaviorCountrySimMK.R')

####################################################################################################
## Simulate Tanzania but with assortativity/heterogeneity & varied amounts of pre-couple extra-couple
####################################################################################################
## cc <- which(ds.nm=='Tanzania') 
countries <- 1:length(ds.nm)
each.val <- 200 ##  equates to ~100,000 couples

blocks <- expand.grid(group = countries,
                      het.beh=T, het.beh.sd=1:3, het.beh.cor=sqrt(c(0,.25,.5,.75)),
                      het.gen=F, het.gen.sd=0, het.gen.cor=0,                                    
                      bmb.sc = exp(seq(log(.1), log(10), l=15)))
blocks <- rbind(blocks, expand.grid(group = countries,
                      het.beh=T, het.beh.sd=0, het.beh.cor=0,
                      het.gen=F, het.gen.sd=0, het.gen.cor=0,                                    
                      bmb.sc = exp(seq(log(.1), log(10), l=15))))
blocks$bfb.sc <- blocks$bme.sc <- blocks$bfe.sc <- blocks$bmb.sc ## all contact coefficients scaled same way
blocks$bmp.sc <- blocks$bfp.sc <- 1
blocks <- blocks[with(blocks, order(group, het.beh.sd, het.beh.cor, bmb.sc)),]

blocks$jobnum <- 1:nrow(blocks)
nn <- nrow(blocks)
head(blocks,50)

outdir <- file.path('results','AssortBehaviorCountrySim')
if(!file.exists(outdir))      dir.create(outdir) # create directory if necessary
batchdirnms <- rep(NA, length(countries))
for(cc in 1:length(countries)) {
    batchdirnms[cc] <- file.path(outdir, ds.nm[cc]) # setup directory for country
    if(!file.exists(batchdirnms[cc]))      dir.create(batchdirnms[cc]) # created if necessary
    if(!file.exists(file.path(batchdirnms[cc],'routs')))      dir.create(file.path(batchdirnms[cc], 'routs')) # setup directory to store Rout files
}

load(file=file.path(outdir, 'JobsToDo.Rdata')) ## for finishing up jobs from last run that didn't get finished due to cluster problems.

num.doing <- 0
sink("CountryAssortSim.txt")         # create a control file to send to the cluster
## LEFT OFF HERE
for(ii in 1:nn) {
    jb <- ii                   # job num
    cmd <- with(blocks, {
        paste("R CMD BATCH '--args jobnum=", ii, " batchdirnm=\"", batchdirnms[group[ii]], "\"", " nc=", nc,
              " group.ind=", group[ii], " substitute=F sub.betas=F counterf.betas=F",
              " s.epic=", group[ii],  " s.demog=", group[ii],
              " s.bmb=", group[ii], " s.bfb=", group[ii],
              " s.bme=", group[ii], " s.bfe=", group[ii],
              " s.bmp=", group[ii], " s.bfp=", group[ii], 
              " death=T",
              " acute.sc=7 late.sc=5 aids.sc=0 dur.ac=2 dur.lt=10 dur.aids=10", # acute phase varying throughout loop
              " bmb.sc=", bmb.sc[ii], " bfb.sc=", bfb.sc[ii],
              " bme.sc=", bme.sc[ii], " bfe.sc=", bfe.sc[ii],
              " bmp.sc=", bmp.sc[ii], " bfp.sc=", bfp.sc[ii],
              " het.b=F het.b.sd=0 het.b.cor=0",
              " het.e=F het.e.sd=0 het.e.cor=0",
              " het.p=F het.p.sd=0 het.p.cor=0",
              " het.gen=", het.gen[ii], " het.gen.sd=", het.gen.sd[ii], " het.gen.cor=", het.gen.cor[ii],
              " het.beh=", het.beh[ii], " het.beh.sd=", het.beh.sd[ii], " het.beh.cor=", het.beh.cor[ii],
              " hilo=F phihi=.2 phi.m=.2 phi.f=.2 rrhi.m=10 rrhi.f=10", ## default hilo parameters, not used
              " scale.by.sd=T scale.adj=1",
              " infl.fac=200 maxN=10^5 sample.tmar=F",
              " psNonPar=F seed=1 tmar=(65*12):(113*12) each=", each.val,
              " tint=113*12' SimulationStarter.R ", file.path(batchdirnms[group[ii]], "routs", paste0(ds.nm[group[ii]], ii, ".Rout")), sep='')
    })
                                        #     if(totn %in% jtd & ii %in% 89:92 & aa==7) { ## for finishing up jobs that didn't get properly submitted (cluster issues sometimes)
    if(blocks$group[ii]==12 & ii %in% jtd) { ## for finishing up jobs that didn't get properly submitted (cluster issues sometimes)
        num.doing <- num.doing+1
        cat(cmd)               # add command
        cat('\n')              # add new line
    }
}
sink()
print(nn)
print(num.doing)
save(blocks, file = file.path(outdir,'blocks.Rdata')) # these are country-acute phase specific blocks
####################################################################################################

