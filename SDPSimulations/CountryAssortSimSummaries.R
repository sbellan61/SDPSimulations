####################################################################################################
## Collect, summarize, & visualize results from simulations of one
## country with varying behavior, assortativity & heterogeneity.
####################################################################################################
rm(list=ls())                           # clear workspace
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02413/sbellan/DHSProject/SDPSimulations/')
source('PlotFunctions.R')                    # load functions to collect & plot results
source('SimulationFunctions.R')                   # load simulation functions
library(abind)                          # array binding
load('data files/ds.nm.all.Rdata')        # country names
load('data files/dframe.s.Rdata')        # SDP by survey
do.again <- T                           # collect results again (otherwise load cfs.Rdata)
show.pts <- T                           # show observed SDP in real DHS data
## source('CountryAssortSimSummaries.R')

dir.results <- file.path('results','AssortBehaviorCountrySim')
dir.figs <- file.path(dir.results, 'Figures')        # make a directory to store figures
if(!file.exists(dir.figs)) dir.create(dir.figs)      # create it
## load 'blocks' which gives info on all simulations within each country-acute group
load(file.path(dir.results, 'blocks.Rdata'))
## Load all results files in the results directory (all Rdata files except blocks & cfs)
fs <- list.files(pattern = '.Rdata', path = file.path(dir.results), recursive = T, full.names = T)
fs <- fs[!grepl('blocks',fs) & !grepl('cfs',fs) &  !grepl('JobsToDo',fs)]
#fs <- fs[grepl('Acute7',fs)]
print(length(fs))

## coll: collects all results into data frame of input parameters, and array of time series
if(!file.exists(file.path(dir.results, 'cfs.Rdata')) | do.again) {
  cfs <- coll(fs, nc = 12, give.ev = F, lbrowse=F, trace=F, browse=F)
  cfs$cframe <- cfs$cframe[order(cfs$cframe$job),] # order by jobs
  cfs$t.arr <- cfs$t.arr[,,order(cfs$cframe$job)]          # ditto
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
nac <- 8 #length(acutes)                   # how many
col.pl <- 'black'                       # base plot color
mend <- max(blocks$end)                 # last job of each block
countries <- unique(cframe$group.ind)   # countries simulated
countries <- countries[order(countries)]
ngroup <- length(countries)             # how many
jtd.all <- with(blocks, jobnum[group==12])
jtd <- jtd.all[!jtd.all %in% cframe$job]
print(paste("didn't do jobs:",paste(head(jtd,50), collapse=','))) # check to see if any jobs didn't complete
save(jtd, file=file.path(dir.results, 'JobsToDo.Rdata'))

####################################################################################################
## Need to fit these simulations with a homogenous model

names(cframe)
blocks
