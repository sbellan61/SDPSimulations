####################################################################################################
## Makes control files for each analysis within which each line giving one R CMD BATCH command line
## to run on a cluster.
####################################################################################################
rm(list=ls())                                  # clear workspace
load("data files/ds.nm.all.Rdata") # country names
load('data files/pars.arr.ac.Rdata')    # load acute phase relative hazards used to fit (in.arr[,,2])
load('data files/CFJobsToDo.Rdata') ## for finishing up jobs from last run that didn't get finished due to cluster problems.
hazs <- c('bmb','bfb','bme','bfe','bmp','bfp') #  transmission coefficient names, for convenience
nc <- 12                                       # core per simulation
## source('AssortMK.R')

####################################################################################################
####################################################################################################
##   COUNTERFACTUAL ANALYSIS
####################################################################################################
## Using an HIV transmission model fit to Demographic and Health Surveys in 18 sub-Saharan African
## countries, we used counter-factual simulations to examine how serodiscordant proportions are
## affected by assortativity between patners
####################################################################################################
countries <- 1:length(ds.nm)
countries <- which(ds.nm=='Zambia')
cc <- which(ds.nm=='Uganda') ## to get Ugandan prameters from DHS
each <- 200 ##  equates to ~100,000 couples
blocks <- expand.grid(group = countries,
                      acute.sc = 7,
                      late.sc = 1, aids.sc = 1, death = T,
                      het.gen.sd = seq(0,3, by = .5),
                      het.gen.cor = c(0,.5,.7, .9),
                      bmb.sc = 1, bfb.sc = 1,
                      bme.sc = 1, bfe.sc = 1,
                      bmp.sc = 1, bfp.sc = 1)
blocks$het.gen <- blocks$het.gen.sd > 0

counterf.betas <- F ## change betas in counterfactuals? if not change beta_within & c's (so beta_within affects all routes)
sub.betas <- F      ## substitute betas? if not beta_within & c's
substitute <- F ## not a substitution analysis
totn <- 0 ## total number of simulations (steps up to final value)
num.doing <- 0
maxN <- 10^5 ## max pseudopopulation size
sample.tmar <- F ## sample marital (couple formation) date from copulas?
psNonPar <- F ##  use non-parametric couple pseudo-population builder?
outdir <- file.path('results','CounterAssort')
if(!file.exists(outdir))      dir.create(outdir) ## create directory if necessary

sink("Assort.txt") ## create a control file to send to the cluster
for(ii in 1:nrow(blocks)) { ## for each simulation
  acdirnm <- file.path(outdir,paste0('Acute', blocks$acute.sc[ii])) ## set up directory for each acute phase relative hazard
  if(!file.exists(acdirnm))      dir.create(acdirnm) ## create directory if necessary
  batchdirnm <- file.path(acdirnm, ds.nm[blocks$group[ii]]) ## setup directory for country
  if(!file.exists(batchdirnm))      dir.create(batchdirnm) ## created if necessary
  if(!file.exists(file.path(batchdirnm,'routs')))      dir.create(file.path(batchdirnm, 'routs')) ## setup directory to store Rout files
  jb <- ii                   # job num
  totn <- totn+1             # total jobs
  cmd <- with(blocks[ii,], 
              paste("R CMD BATCH '--args jobnum=", totn, " simj=", ii, " batchdirnm=\"", batchdirnm, "\"", " nc=", nc,
                    " group.ind=", group, " substitute=", substitute, " sub.betas=", sub.betas, " counterf.betas=", counterf.betas,
                    " death=", death, " acute.sc=", acute.sc, " late.sc=", late.sc," aids.sc=", aids.sc, 
                    " bmb.sc=", bmb.sc, " bfb.sc=", bfb.sc,
                    " bme.sc=", bme.sc, " bfe.sc=", bfe.sc,
                    " bmp.sc=", bmp.sc, " bfp.sc=", bfp.sc,
                    " s.bmb=", group, " s.bfb=", group,
                    " s.bme=", group, " s.bfe=", group,
                    " s.bmp=", group, " s.bfp=", group,
                    " s.demog=", group, " s.epic=", group,
                    " het.gen=", het.gen, " het.gen.sd=", het.gen.sd, " het.gen.cor=", het.gen.cor,
                    " sample.tmar=", sample.tmar,
                    " seed=1 tmar=(65*12):(113*12) each=", each, " maxN=", maxN, 
                    " tint=113*12' SimulationStarter.R ", file.path(batchdirnm, "routs", paste0(ds.nm[group], ii, ".Rout")), sep='')
              )
  num.doing <- num.doing+1
  cat(cmd)               # add command
  cat('\n')              # add new line
}
sink()
save(blocks, file = file.path(outdir,'blocks.Rdata')) # these are country-acute phase specific blocks
## ##################################################################################################
print(totn)
print(num.doing)
