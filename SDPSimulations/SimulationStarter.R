######################################################################
## Drivers of serodiscordance patterns and sub-Saharan Africa
###################################################################### 
## Steve Bellan, 2013-2017
## steve.bellan@gmail.com
#################################################################################################### 
##  The goal of this analysis is to assess the extent to which pre-, extra-, & within-couple
##  transmission, AIDS mortality, and acute phase infectivity drive observed levels of
##  serodiscordance in sub-Saharan Africa as well as to determine which country-specific
##  characteristics can explain their variable serodiscordance proportions and epidemic intensities.
#################################################################################################### 
## This code is called by R CMD BATCH with input arguments (each line of a control file sent to a
## cluster), runs a simulation using functions from sim fxns3.R, and then saves the output in a
## specified directory structure
#################################################################################################### 
rm(list=ls())                           # clear workspace

## acute.sc=1;country=15;batch=1;simj=1;lab="as;fitted";death="TRUE";het.gen="FALSE";het.gen.sd=0;het.gen.cor=0;het.beh="FALSE";het.beh.sd=0;het.beh.cor=0;het.b="FALSE";het.b.sd=0;het.b.cor=0;het.e="FALSE";het.e.sd=0;het.e.cor=0;het.p="FALSE";het.p.sd=0;het.p.cor=0;bmb.sc=1;bfb.sc=1;bme.sc=1;bfe.sc=1;bmp.sc=1;bfp.sc=1;late.sc=1;aids.sc=1;group=15;s.epic=15;s.demog=15;scale.by.sd="TRUE";scale.adj=1;infl.fac=200;maxN=1e+05;sample.tmar="FALSE";psNonPar="FALSE";each=200;jobnum=57;seed=1;out.dir="results/CounterFactual";sim.nm="CF";substitute=F;tmar="(65*12):(113*12)";tint=113*12

args=(commandArgs(TRUE))                # load arguments from R CMD BATCH 
if(length(args)>0)  {## Then cycle through each element of the list and evaluate the expressions.
    for(i in 1:length(args)) {
        eval(parse(text=args[[i]]))
    }  }
tmar <- eval(parse(text=tmar)) ## convert string to series
set.seed(seed)
vfreq <- 200 ## how often to report on results
source("SimulationFunctions.R")                   # load simulation functions from script
country <- group              # set country to country-group index
odat <- dat          # backup original data in this workspace before it gets manipulated
print(ds.nm[country])                   # Print country name.

######################################################################
## Transmission coefficients
## Load parameters for acute phase assuming (acute.sc)
pars.arr <- out.arr[,,which(in.arr[,1,2]==acute.sc),] # median and credible intervals for transmission coefficient estimates for this acute relative hazard
hazs <- c("bmb","bfb","bme","bfe","bmp","bfp") # six gender-route specific transmission coefficients (*b*efore-, *e*xtra-, from-*p*artner-) for *m*ale & *f*emale
spars <- pars.arr[hazs,2,country]              # get transmission coefficients from base country
## If substituting, substitute parameters/epidemic curves out for those from donor country
if(substitute) {
    ##  For each country substitute things from others and see how close
    ##  it gets to the serodiscordance levels of the other country.
######################################################################
    ##  Epidemic curve
    if(s.epic!=country) {     # if substituting, replce it in data (otherwise it defaults to original)
        s.epic.nm <- ds.nm[s.epic]              # epidemic curve to use (country name)
        s.epic.ind <- which(colnames(epicf)==s.epic.nm) # epidemic curve to use (column index of epicf/m matrices)
        ## substitute epidemic curves into data
        dat$epic.ind <- s.epic.ind
        dat$epic.nm <- s.epic.nm
    }
    if(sub.betas) { ## substituting betas between countries
        ##  substitute parameters from substitution (donor) country. At most two of these substitutions will change spars.
        for(hh in hazs[1:4]) {
            s.hh <- get(paste('s.',hh))
            spars[hh] <- pars.arr[hh,2,s.hh]
        }
        ## substitute demography is just done by using parameters & epi curves from the opposite country in psrun() below.
    }else{ ## substituting HIV transmission rate & contact coefficients
        for(hh in hazs[1:4]) {
            s.hh <- get(paste('s.',hh))
            chh <- sub('b', 'c',hh)
            spars[hh] <- pars.arr[hh,2,s.hh] * pars.arr[chh,2,s.hh]
        }
    }
    for(hh in hazs[4:6]) { ##w/in couple transmission doesn't vary based on sub.betas
        s.hh <- get(paste('s.',hh))
        spars[hh] <- pars.arr[hh,2,s.hh]
    }
}else{
    s.epic.nm <- NA # not substituting epidemic curves, this will cause rcop() to use default country epidemic curves
    s.epic.ind <- NA
}

parmArgs <- subsArgs(as.list(environment()), psrun)
parmArgs$pars <- spars
## parmArgs$s.epic.nm <- s.epic.nm ## not used below for some reason?
## parmArgs$s.epic.ind <- s.epic.ind
parmArgs$each <- 1

## Simulate couples transmission model (calling psrun() from sim fxns3.R). Output (temp) is the name of the file that is produced.
temp <- do.call(psrun, args=parmArgs)

load(temp)                              # load output of simulations
names(output)                           # list objects summarized below
## jobnum=which job; evout=line list of each couple; tss=time series of pseudopopulation, pars=input
## parameters; tmar=marriage cohort dates; each=# couples per marriage cohort
## 
head(output$evout)                      # columns explained below
## with(output$evout, cor(log(m.het.gen), log(f.het.gen)))
## with(output$evout, range(log(m.het.gen), log(f.het.gen)))
## xtabs(~m.het.hilo + f.het.hilo, output$evout) / nrow(output$evout)
## xtabs(~m.het.hilo + f.het.hilo + ser, output$evout)

## uid=unique couple identifier; ds=dat set; ser=couple serostatus (1:4::{++,+-,-+,--});
## tms/tfs=male/female sexual debuts (months since 1900); tmar=couple formation (i.e. marriage)
## date; tint=DHS interview/testing date; mardur.mon=couple duration in months; circ=circumciscion
## status of male; mage/fage=male/female age at interview; epic.ind=epidemic curve index in
## allepicm/f; epic.nm=epidemic curve name; group=analysis group; m/fser=m/f serostatus at
## interview; f/mdoi=m/f date of infection; m/fdod=m/f date of death; m/fcoi=m/f cause (i.e. route)
## of infection; m/f.het.gen/beh/b/e/p=lognormal risk deviate for m/f for genetic, behavioral or
## route-specific heterogeneity; m/fcoi.phase=partner's phase (acute/chronic/late) for infection;
## tend=end of couple-time (interview or death); dm/fage=date at which m/f age out of DHS cohort
## (50yr for f, 60yr for m); taend=end of couple-time (interview, aging out of DHS cohort, or
## death); m/falive=m/f alive at interview; tmsdc=time spent as M+F- couple; tfsdc=time spent as
## M-F+ couple; tccc=time spent ++; m/fsdc=ever M/F+ serodicsordant?; ccc=ever ++?
tail(output$tss)                        # columns explained below
## tt=months since 1900; yr=year; b.ss=#couples -- *b*efore couple formation; b.mm/b.ff=# couples
## +-/-+ before cf; b.hh=#couples ++ before cf; ss/ff/mm/hh=# couples in each serostatus formed &
## still alive; d.xx=# couples in each serostatus but >=1 partner dead; a.xx=# couples in each
## serostatus but >=1 partner aged out of DHS cohort (gender that's aged is appended); alive=#
## couples alive; inf.alive=# couples infected & alive; inf=# couples infected; cu.xx=#cumulative
## infections by gender-route combination
print(output$pars)                      # input parameters
print(spars)                            # inputted hazards
