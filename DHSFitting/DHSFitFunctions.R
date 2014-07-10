## Pre-couple transmission iterator
pre.couple <- function(seros.active, pmb, pfb, pmb.a, pfb.a, uncond.mort=T) {
    if(class(seros.active)=='numeric')  seros.active <- t(as.matrix(seros.active))
    tp <- seros.active ## old temporary array to update from
    seros.active[,'s..']   <- tp[,'s..'] * (1-pmb) * (1-pfb)
    ## transmission and alive, these states are used for fitting
    seros.active[,'mb.a1A'] <- tp[,'s..'] * pmb.a * (1-pfb)
    seros.active[,'mb.a2A'] <- tp[,'mb.a1A']*(1-pfb)
    seros.active[,'mb.A'] <- tp[,'mb.a2A']*(1-pfb) + tp[,'mb.A']*(1 - pfb)
    seros.active[,'f.ba1A'] <- tp[,'s..']* pfb.a *(1-pmb)
    seros.active[,'f.ba2A'] <- tp[,'f.ba1A']*(1 - pmb)
    seros.active[,'f.bA'] <- tp[,'f.ba2A']*(1 - pmb) + tp[,'f.bA']*(1 - pmb)
    p.mfirst.a <- pmb.a / (pmb.a+pfb.a)
    p.ffirst.a <- 1-p.mfirst.a
    p.mfirst.a[is.na(p.mfirst.a)] <- 0
    p.ffirst.a[is.na(p.ffirst.a)] <- 0                
    seros.active[,'hb1b2A'] <- tp[,'hb1b2A'] + p.mfirst.a * tp[,'s..']* pmb.a * pfb.a + (tp[,'mb.a1A'] + tp[,'mb.a2A'] + tp[,'mb.A']) * pfb.a
    seros.active[,'hb2b1A'] <- tp[,'hb2b1A'] + p.ffirst.a * tp[,'s..']* pmb.a * pfb.a + (tp[,'f.ba1A'] + tp[,'f.ba2A'] + tp[,'f.bA']) * pmb.a
    if(uncond.mort) { ## if calculating probabilities of transmission with or without death (not used for fitting)
        seros.active[,'mb.a1'] <- tp[,'s..'] * pmb * (1-pfb)
        seros.active[,'mb.a2'] <- tp[,'mb.a1'] * (1 - pfb)
        seros.active[,'mb.'] <- tp[,'mb.a2'] * (1 - pfb) + tp[,'mb.'] * (1 - pfb)
        seros.active[,'f.ba1'] <- tp[,'s..'] * pfb * (1-pmb)
        seros.active[,'f.ba2'] <- tp[,'f.ba1'] * (1 - pmb)
        seros.active[,'f.b'] <- tp[,'f.ba2'] * (1 - pmb) + tp[,'f.b'] * (1 - pmb)
        p.mfirst <- pmb / (pmb+pfb)
        p.ffirst <- 1-p.mfirst
        p.mfirst[is.na(p.mfirst)] <- 0
        p.ffirst[is.na(p.ffirst)] <- 0                
        seros.active[,'hb1b2'] <- tp[,'hb1b2'] + p.mfirst  *  tp[,'s..'] * pmb * pfb + (tp[,'mb.a1'] + tp[,'mb.a2'] + tp[,'mb.'])  *  pfb
        seros.active[,'hb2b1'] <- tp[,'hb2b1'] + p.ffirst  *  tp[,'s..'] * pmb * pfb + (tp[,'f.ba1'] + tp[,'f.ba2'] + tp[,'f.b'])  *  pmb
    }
    return(seros.active)
}

## Within-couple transmission iterator
within.couple <- function(seros.active, pme, pfe, pmp, pfp, pmp.ac, pfp.ac, pme.a, pfe.a, pmp.a, pfp.a, pmp.a.ac, pfp.a.ac, uncond.mort=T) {
    if(class(seros.active)=='numeric')  seros.active <- t(as.matrix(seros.active))
    tp <- seros.active
    seros.active[,'s..'] <- tp[,'s..']*(1-pme)*(1-pfe)
    ## ##################################################
    ## Joint with survival
    seros.active[,'mb.a1A'] <- 0
    seros.active[,'mb.a2A'] <- tp[,'mb.a1A']*(1-pfe)*(1-pfp.ac)
    seros.active[,'mb.A']  <- tp[,'mb.a2A']*(1-pfe)*(1-pfp.ac) + tp[,'mb.A']*(1-pfe)*(1-pfp) 
    seros.active[,'me.a1A'] <- tp[,'s..']*pme.a*(1-pfe)
    seros.active[,'me.a2A'] <- tp[,'me.a1A']*(1-pfe)*(1-pfp.ac)
    seros.active[,'me.A']  <- tp[,'me.a2A']*(1-pfe)*(1-pfp.ac) + tp[,'me.A']*(1-pfe)*(1-pfp)
    seros.active[,'f.ba1A'] <- 0
    seros.active[,'f.ba2A'] <- tp[,'f.ba1A']*(1-pme)*(1-pmp.ac)
    seros.active[,'f.bA']  <- tp[,'f.ba2A']*(1-pme)*(1-pmp.ac) + tp[,'f.bA']*(1-pme)*(1-pmp)
    seros.active[,'f.ea1A'] <- tp[,'s..']*pfe.a*(1-pme)
    seros.active[,'f.ea2A'] <- tp[,'f.ea1A']*(1-pme)*(1-pmp.ac)
    seros.active[,'f.eA']  <- tp[,'f.ea2A']*(1-pme)*(1-pmp.ac) + tp[,'f.eA']*(1-pme)*(1-pmp)
    ##              hb1b2A'] <- hb1b2A'] # Doesn't change during couple duration
    ##              hb2b1A'] <- hb2b1A'] # Doesn't change during couple duration                
    seros.active[,'hbeA'] <- tp[,'hbeA']  + (tp[,'mb.a1A'] + tp[,'mb.a2A'])*(1-pfp.ac)*pfe.a + tp[,'mb.A']*(1-pfp)*pfe.a
    seros.active[,'hebA'] <- tp[,'hebA']  + (tp[,'f.ba1A'] + tp[,'f.ba2A'])*(1-pmp.ac)*pme.a + tp[,'f.bA']*(1-pmp)*pme.a
    seros.active[,'hbpaA'] <- tp[,'hbpaA'] + (tp[,'mb.a1A'] + tp[,'mb.a2A'])*pfp.a.ac
    seros.active[,'hbpA'] <- tp[,'hbpA']  + tp[,'mb.A']*pfp.a                
    seros.active[,'hpbaA'] <- tp[,'hpbaA'] + (tp[,'f.ba1A'] + tp[,'f.ba2A'])*pmp.a.ac
    seros.active[,'hpbA'] <- tp[,'hpbA']  + tp[,'f.bA']*pmp.a                
    seros.active[,'hepaA'] <- tp[,'hepaA'] + (tp[,'me.a1A'] + tp[,'me.a2A'])*pfp.a.ac
    seros.active[,'hepA'] <- tp[,'hepA']  + tp[,'me.A']*pfp.a                
    seros.active[,'hpeaA'] <- tp[,'hpeaA'] + (tp[,'f.ea1A'] + tp[,'f.ea2A'])*pmp.a.ac
    seros.active[,'hpeA'] <- tp[,'hpeA']  + tp[,'f.eA']*pmp.a                
    ## for individuals infected in the same month, assign
    ## the order of infection based on competing risks
    ## formula, but if the denominator is 0, replace both
    ## with 0 to avoid errors.
    p.mfirst.a <- pme.a / (pme.a+pfe.a)
    p.ffirst.a <- 1-p.mfirst.a
    p.mfirst.a[is.na(p.mfirst.a)] <- 0
    p.ffirst.a[is.na(p.ffirst.a)] <- 0                
    seros.active[,'he1e2A'] <- tp[,'he1e2A'] + p.mfirst.a * tp[,'s..']*pme.a*pfe.a + (tp[,'me.a1A'] + tp[,'me.a2A'])*(1-pfp.ac)*pfe.a + tp[,'me.A']*(1-pfp)*pfe.a
    seros.active[,'he2e1A'] <- tp[,'he2e1A'] + p.ffirst.a * tp[,'s..']*pme.a*pfe.a + (tp[,'f.ea1A'] + tp[,'f.ea2A'])*(1-pmp.ac)*pme.a + tp[,'f.eA']*(1-pmp)*pme.a
    if(uncond.mort) { ## if calculating probabilities of transmission with or without death (not used for fitting)
        seros.active[,'mb.a1'] <- 0 
        seros.active[,'mb.a2'] <- tp[,'mb.a1'] * (1-pfe)*(1-pfp.ac)
        seros.active[,'mb.']   <- tp[,'mb.a2'] * (1-pfe)*(1-pfp.ac) + tp[,'mb.'] * (1-pfe)*(1-pfp)            
        seros.active[,'me.a1'] <- tp[,'s..'] * pme*(1-pfe)
        seros.active[,'me.a2'] <- tp[,'me.a1']*(1-pfe)*(1-pfp.ac) 
        seros.active[,'me.']   <- tp[,'me.a2']*(1-pfe)*(1-pfp.ac) + tp[,'me.']*(1-pfe)*(1-pfp) 
        seros.active[,'f.ba1'] <- 0
        seros.active[,'f.ba2'] <- tp[,'f.ba1']*(1-pme)*(1-pmp.ac)
        seros.active[,'f.b']   <- tp[,'f.ba2']*(1-pme)*(1-pmp.ac) + tp[,'f.b']*(1-pme)*(1-pmp)            
        seros.active[,'f.ea1'] <- tp[,'s..']*pfe*(1-pme)
        seros.active[,'f.ea2'] <- tp[,'f.ea1']*(1-pme)*(1-pmp.ac)
        seros.active[,'f.e']   <- tp[,'f.ea2']*(1-pme)*(1-pmp.ac) + tp[,'f.e']*(1-pme)*(1-pmp) ## hb1b2/b2b1 not here b/c don't change during marriage
        seros.active[,'hbe']  <- tp[,'hbe'] + (tp[,'mb.a1']+ tp[,'mb.a2'])*(1-pfp.ac)*pfe + tp[,'mb.']*(1-pfp)*pfe
        seros.active[,'heb']  <- tp[,'heb'] + (tp[,'f.ba1']+ tp[,'f.ba2'])*(1-pmp.ac)*pme + tp[,'f.b']*(1-pmp)*pme
        seros.active[,'hbpa'] <- tp[,'hbpa']+ (tp[,'mb.a1']+ tp[,'mb.a2'])*pfp.ac
        seros.active[,'hbp']  <- tp[,'hbp'] + tp[,'mb.']*pfp            
        seros.active[,'hpba'] <- tp[,'hpba']+ (tp[,'f.ba1']+ tp[,'f.ba2'])*pmp.ac
        seros.active[,'hpb']  <- tp[,'hpb'] + tp[,'f.b']*pmp
        seros.active[,'hepa'] <- tp[,'hepa']+ (tp[,'me.a1']+ tp[,'me.a2'])*pfp.ac
        seros.active[,'hep']  <- tp[,'hep'] + tp[,'me.']*pfp            
        seros.active[,'hpea'] <- tp[,'hpea']+ (tp[,'f.ea1']+ tp[,'f.ea2'])*pmp.ac
        seros.active[,'hpe']  <- tp[,'hpe'] + tp[,'f.e']*pmp            
        p.mfirst <- pme / (pme+pfe)
        p.ffirst <- 1-p.mfirst
        p.mfirst[is.na(p.mfirst)] <- 0
        p.ffirst[is.na(p.ffirst)] <- 0                
        seros.active[,'he1e2'] <- tp[,'he1e2'] + p.mfirst * tp[,'s..']*pme*pfe + (tp[,'me.a1'] + tp[,'me.a2'])*(1-pfp.ac)*pfe + tp[,'me.']*(1-pfp)*pfe
        seros.active[,'he2e1'] <- tp[,'he2e1'] + p.ffirst * tp[,'s..']*pme*pfe + (tp[,'f.ea1'] + tp[,'f.ea2'])*(1-pmp.ac)*pme + tp[,'f.e']*(1-pmp)*pme
    }
    return(seros.active)
}


## Get (1) logical matrices for whether male, female, or either partner is sexually active for each
## tt-th month of pre-couple transmission, (2) country-level prevalence for each couple's tt-th month of before-couple duration, (3) probability either partner would survive to interview date if infected in tt-th month of pre-couple transmission.
pre.prep <- function(dat) {
    out.names <- c('m.sex','f.sex','e.sex','pre.fprev','pre.mprev','pre.msurv','pre.fsurv')
    for(ob in out.names[1:3]) assign(ob, matrix(NA, nrow(dat), max(dat$bd)))
    for(ob in out.names[4:length(out.names)]) assign(ob, matrix(0, nrow(dat), max(dat$bd))) ## no transmission occurs by default when they aren't sexually active
    for(tt in 1:max(dat$bd)) { ## for each month in the before-couple duration (bd)
        m.sex[,tt] <- dat$tmar-dat$bd+tt-1 >= dat$tms & dat$tmar-dat$bd+tt-1 < dat$tmar
        f.sex[,tt] <- dat$tmar-dat$bd+tt-1 >= dat$tfs & dat$tmar-dat$bd+tt-1 < dat$tmar
        e.sex[,tt] <- m.sex[,tt] |f.sex[,tt] ## either are active 
        ## opposite gender population prevalence in month tt of bd
        pre.fprev[m.sex[,tt],tt] <- epicf[cbind(dat$tmar[m.sex[,tt]]-dat$bd[m.sex[,tt]]+tt-1, dat$epic.ind[m.sex[,tt]])]
        pre.mprev[f.sex[,tt],tt] <- epicm[cbind(dat$tmar[f.sex[,tt]]-dat$bd[f.sex[,tt]]+tt-1, dat$epic.ind[f.sex[,tt]])]
        ## probability each partner survives to interview date given they are infected in month tt of bd
        pre.msurv[m.sex[,tt],tt] <- csurv[cbind(dat$mage[m.sex[,tt]]-dat$cd[m.sex[,tt]]-dat$bd[m.sex[,tt]]+tt-1, dat$cd[m.sex[,tt]]+dat$bd[m.sex[,tt]]-tt+1)]
        pre.fsurv[f.sex[,tt],tt] <- csurv[cbind(dat$fage[f.sex[,tt]]-dat$cd[f.sex[,tt]]-dat$bd[f.sex[,tt]]+tt-1, dat$cd[f.sex[,tt]]+dat$bd[f.sex[,tt]]-tt+1)]
    }
    pre.prepout <- list(m.sex, f.sex, e.sex, pre.fprev, pre.mprev, pre.msurv, pre.fsurv)
    names(pre.prepout) <- out.names
    return(pre.prepout)
}
## Get (1) logical matrices for whether couple is active for each tt-th month of marital durations
## (i.e. how long was the couple duration), (2) country-level prevalence for each couple's tt-th
## month of couple duration, (3) probability either partner would survive to interview date if
## infected in tt-th month of within-couple transmission, (4) ART coverage for tt-the month of
## couple duration if reducing within-couple transmission by ART coverage.
within.prep <- function(dat) {
    out.names <- c('fmd','within.fprev','within.mprev', 'within.art.cov', 'within.msurv', 'within.fsurv')
    for(ob in out.names) assign(ob, matrix(NA, nrow(dat), max(dat$cd)))
    for(tt in 1:max(dat$cd-1)) { ## Loop through marital duration from 1st month to last month of longest marriage
        ff <- dat$cd >= tt ## are partners formed in a couple?
        within.fprev[ff,tt] <- epicf[cbind(dat$tmar[ff]+(tt-1), dat$epic.ind[ff])]
        within.mprev[ff,tt] <- epicm[cbind(dat$tmar[ff]+(tt-1), dat$epic.ind[ff])]
        within.art.cov[ff,tt] <- art.cov[cbind(dat$tmar[ff]+(tt-1), dat$epic.ind[ff])]
        within.msurv[ff,tt] <- csurv[cbind(dat$mage[ff]-dat$cd[ff]+tt-1, dat$tint[ff] - (dat$tmar[ff] + tt - 1))]
        within.fsurv[ff,tt] <- csurv[cbind(dat$fage[ff]-dat$cd[ff]+tt-1, dat$tint[ff] - (dat$tmar[ff] + tt - 1))]
        fmd[,tt] <- ff        
    }
    within.prepout <- list(fmd,within.fprev,within.mprev, within.art.cov, within.msurv, within.fsurv)
    names(within.prepout) <- out.names
    return(within.prepout)
}

## Calculate likelihood
pcalc <- function(pars, dat, browse = F, 
                  give.ser = F, ## return serostatus probability matrix
                  acute.sc = 7, ## acute phase infectiousness
                  ## using global variables for these
                  ## partner.arv = F, ## ART coverage affects within-partnership transmission?
                  ## low.coverage.arv = F, ## if T, only 50% of those on ART are not infectious
                  uncond.mort = T, ## track serostate probabilities unconditional on mortality (not currently working)
                  survive = T, ## include mortality in the model (otherwise fitting to death-unconditional probability states)
                  lrho.sd = 1/2) ## sd on lrho prior
    {
        cov.scalar <- ifelse(low.coverage.arv, .5, 1) ## if assuming only 50% of those on ART are not-infectious
        K <- nrow(dat)
        ser.nms <- c('hh','mm','ff','ss')
        if(sum(pars[1:5]<0)>0) { ## if any parameters are <0 
            lprob <- -Inf ## then the model must be rejected so we return logprob =-Inf
            empty.vars <- c('probs','pop.avs','proj12','pser.a','pser','rrs')
            for(ii in 1:length(empty.vars)) assign(empty.vars[ii], NA) ## set other outputs to NA
        }else{
            for(ii in 1:length(pars)) assign(names(pars)[ii], as.numeric(pars[ii])) ## get pars from pars: bmb, bfb, bme, bfe, bmp, lrho
            rho <- exp(lrho) # feeding in log(rho) log(m->f / f->m) transmission rates
            bfp <- bmp * rho
            ## Assign state variable names
            state.var.nms.uncond <- c('s..',
                                      'mb.a1', 'mb.a2', 'mb.', 'me.a1', 'me.a2', 'me.', #M+F- by routes & acute phase (a)
                                      'f.ba1', 'f.ba2', 'f.b', 'f.ea1', 'f.ea2', 'f.e', #M-F+ by routes & acute phase (a)
                                      'hb1b2', 'hb2b1', 'hbe', 'heb', ## rest are all M+F+
                                      'hbpa', 'hpba', 'hepa', 'hpea', ## acute infected by partner
                                      'hbp', 'hpb', 'hep', 'hpe', ## chronic infected by partner
                                      'he2e1', 'he1e2') ## both extra-couply infected, different orders
            state.var.nms <- c(state.var.nms.uncond, paste0(state.var.nms.uncond[-1],'A')) ## joint with alive at the end
            if(browse) browser()
            seros <- matrix(0, K, length(state.var.nms), dimnames = list(NULL,state.var.nms))
            if(!uncond.mort) seros[,state.var.nms.uncond[-1]] <- NA ## 
            seros[,'s..'] <- 1
            for(tt in 1:max(dat$bd)) { ## for each month in the before-couple duration (bd)
                ## probabilities are non-zero only for times after started having sex and before couple formation
                active <- e.sex[,tt]
                m.haz <- bmb * pre.fprev[active,tt] ## hazards to sexually active men
                f.haz <- bfb * pre.mprev[active,tt] ## hazards to sexually active women
                pmb <- 1 - exp(-m.haz)              ## transmission probabilities
                pfb <- 1 - exp(-f.haz)
                pmb.a <- pmb * pre.msurv[active,tt] ## joint transmission & survival probabilities 
                pfb.a <- pfb * pre.fsurv[active,tt]
                seros[active,] <- pre.couple(seros[active,], pmb=pmb, pfb=pfb, pmb.a=pmb.a, pfb.a=pfb.a, uncond.mort=uncond.mort) ## update serostates
            }
            ## probability of being infected by partner (constant, used inside loop)
            pmp <- 1 - exp(-bmp)
            pfp <- 1 - exp(-bfp)
            pmp.ac <- 1 - exp(-acute.sc * bmp)        ##  during acute
            pfp.ac <- 1 - exp(-acute.sc * bfp)
            ## ##################################################
            for(tt in 1:max(dat$cd-1)) {## Now loop through marriage
                ff <- fmd[,tt]
                m.haz <- bme * within.fprev[ff,tt] ## these vectors are length of active couples (sum(fmd))
                f.haz <- bfe * within.mprev[ff,tt]
                pme <- 1 - exp(-m.haz)
                pfe <- 1 - exp(-f.haz)
                ## adjust probability of being infected by partner by probability partner is infectious (i.e. not on ART)
                if(partner.arv) {
                    within.art.scalar <- (1 - cov.scalar * within.art.cov[ff,tt])
                    pmp <- 1 - exp(-bmp * within.art.scalar)
                    pfp <- 1 - exp(-bfp * within.art.scalar)
                    pmp.ac <- 1 - exp(-acute.sc * bmp * within.art.scalar) ## acute
                    pfp.ac <- 1 - exp(-acute.sc * bfp * within.art.scalar)
                }
                ## Survival probabilities
                s.p.m <- within.msurv[ff,tt]
                s.p.f <- within.fsurv[ff,tt]
                pmp.a <- pmp * s.p.m                ## Transmission probabilities from partner (jointly with survival)
                pfp.a <- pfp * s.p.f
                pmp.a.ac <- pmp.ac * s.p.m                ## acute from partner (jointly with survival)
                pfp.a.ac <- pfp.ac * s.p.f
                pme.a <- pme * s.p.m                ## extra-couple
                pfe.a <- pfe * s.p.f
                ## Update serostates for active couples
                seros[ff,] <- within.couple(seros[ff,], pme=pme, pfe=pfe, pmp=pmp, pfp=pfp, pmp.ac=pmp.ac, pfp.ac = pfp.ac, uncond.mort=uncond.mort,
                                            pme.a = pme.a, pfe.a = pfe.a, pmp.a = pmp.a, pfp.a = pfp.a, pmp.a.ac = pmp.a.ac, pfp.a.ac = pfp.a.ac)
            }
            sum.nms <- c('ss','mm','ff','hh','mmA','ffA','hhA')
            sero.sums <- matrix(NA, K, length(sum.nms), dimnames = list(NULL,sum.nms))
            sero.sums[,'ss'] <- seros[,'s..']
            sero.sums[,'mm'] <- rowSums(seros[,c('mb.a1', 'me.a1', 'mb.a2', 'me.a2', 'mb.', 'me.')])
            sero.sums[,'ff'] <- rowSums(seros[,c('f.ba1', 'f.ea1', 'f.ba2', 'f.ea2', 'f.b', 'f.e')])
            sero.sums[,'hh'] <- rowSums(seros[,c('hb1b2', 'hb2b1', 'hbe', 'heb', 'hbpa', 'hpba', 'hepa', 'hpea', 'hbp', 'hpb', 'hep', 'hpe', 'he1e2', 'he2e1')])
            sero.sums[,'mmA'] <- rowSums(seros[,c('mb.a1A', 'me.a1A', 'mb.a2A', 'me.a2A', 'mb.A', 'me.A')])
            sero.sums[,'ffA'] <- rowSums(seros[,c('f.ba1A', 'f.ea1A', 'f.ba2A', 'f.ea2A', 'f.bA', 'f.eA')])
            sero.sums[,'hhA'] <- rowSums(seros[,c('hb1b2A', 'hb2b1A', 'hbeA', 'hebA', 'hbpaA', 'hpbaA', 'hepaA', 'hpeaA', 'hbpA', 'hpbA', 'hepA', 'hpeA', 'he1e2A', 'he2e1A')])
            if(survive) { ## must NORMALIZE probabilities to 1 for likelihood!
                pser <- sero.sums[,c('hhA', 'mmA', 'ffA', 'ss')]
                probs <- pser[cbind(1:K,dat$ser)] / rowSums(pser) # accounting for survival
            }else{
                pser <- sero.sums[,c('hh', 'mm', 'ff', 'ss')]
                probs <- pser[cbind(1:K,dat$ser)] ## if not accounting for survival
            }
            if(sum(probs==0)==0) { # if non of the serotatuses occur with 0 probability in the current model
                lprob <- sum(log(probs)) + dnorm(log(rho), log(trans.ratio), lrho.sd, log = T) ## add llikelihood and lprior
            }else{ # if some of the serostatuses are 0, then the current parameters have 0 probability
                lprob <- -Inf
            }
        }
        return(list(lprob=lprob,pser=pser,seros=cbind(seros,sero.sums)))
    }

init.fxn <- function(seed = 1)
  {
    set.seed(seed)
    ## sample uniform dispersed on log scale
    lpars <- c(bmb = runif(1, -6, -3),
               bfb = runif(1, -6, -3),
               bme = runif(1, -6, -3),
               bfe = runif(1, -6, -3),
               bmp = runif(1, -6, -3))
    pars <- c(exp(lpars), lrho = runif(1, -1, 1))
    return(pars)
  }

## MCMC SAMPLER
sampler <- function(sd.props = sd.props, inits, dat,
                    acute.sc,
                    multiv = F, covar = NULL, # if multiv, sample from multivariate distribution (calculated during adaptive phase)
                    verbose = T, tell = 100, seed = 1, lrho.sd = 1/2,
                    niter = 6*1000, survive = T, uncond.mort = F,
                    keep.seros = F, ## trace all serostate probabilities
                    nthin = 1,
                    nburn = 1000, browse=F)
  {
    if(browse)  browser()
    set.seed(seed)
    pars <- inits
    vv <- 2
    accept <- 0 ## track each parameters acceptance individually
    cur <- pcalc(pars, acute.sc = acute.sc, dat = dat, survive = survive, uncond.mort = uncond.mort, lrho.sd = lrho.sd)     #calculate first log probability
    lprob.cur <- cur$lprob
    out <- t(as.matrix(c(pars, bfp = as.numeric(pars["bmp"]*exp(pars["lrho"])))))
    if(keep.seros)      seros <- cur$seros else seros <- NULL
    last.it <- 0
    start <- Sys.time()
    while(vv < niter + 1) {
        if(verbose & vv%%tell+1==1) print(paste("on iteration",vv,"of",last.it + niter + 1))
        pars.prop <- pars              #initialize proposal parameterr vector
        ## propose new parameter vector
        if(multiv)          {
            pars.prop <- pars.prop + rmnorm(1, mean = 0, varcov = covar)
            pars.prop <- as.vector(pars.prop) #otherwise is a matrix
            names(pars.prop) <- parnames
          }else{
            pars.prop <- pars.prop + rnorm(length(pars), mean = 0, sd = sd.props)
          }
        ## trace = T if in non-thinned iteration, or the previous one (in case of rejection)
        ## calculate proposal par log probability
        prop <- pcalc(pars.prop, acute.sc = acute.sc, dat = dat, survive = survive)
        lprob.prop <- prop$lprob
        lmh <- lprob.prop - lprob.cur       # log Metropolis-Hastings ratio
        ## if MHR >= 1 or a uniform random # in [0,1] is <= MHR, accept otherwise reject
        if(lmh >= 0 | runif(1,0,1) <= exp(lmh)) {
            pars <- pars.prop
            if(vv>nburn) accept <- accept + 1 #only track acceptance after burn-in
            lprob.cur <- lprob.prop
            cur <- prop
          }
        if(vv%%nthin + 1 ==1) {
            out <- rbind(out,t(as.matrix(c(pars, bfp = as.numeric(pars["bmp"]*exp(pars["lrho"]))))))
            if(keep.seros)      seros <- abind(seros, cur$seros, along = 3)
        }
        vv <- vv+1
    }
    if(verbose) print(paste("took", difftime(Sys.time(),start, units = "mins"),"mins"))
    aratio <- accept/((vv-nburn))
    give <- 1:nrow(out)>(nburn+1)/nthin
    return(list(out = out[give,], aratio = aratio, inits = inits, seros=seros[,,give]))
}

get.rrs <- function(pars.out) {
    with(as.data.frame(pars.out), {
        ## relative rate of transmission coefficient extracouply vs before relationship
        rr.m.eb <- bme/bmb
        rr.f.eb <- bfe/bfb
        rr.m.pe <- bmp/bme
        rr.f.pe <- bfp/bfe
        rr.m.pb <- bmp/bmb        #partner to before
        rr.f.pb <- bfp/bfb
        ## relative rate of transmission coefficient extracouply and before relationship between males and females
        rr.mf.bef <- bmb/bfb
        rr.mf.exc <- bme/bfe
        ## rho is the last one
        ## relative rate of contact/risk paramter (i.e. accounting
        ## for difference in per coital act probability as estimated
        ## from within partnership transmission.
        rr.mf.bef.cont <- rr.mf.bef * exp(lrho)
        rr.mf.exc.cont <- rr.mf.exc * exp(lrho)
        rrs <- data.frame(rr.m.eb = rr.m.eb, rr.f.eb = rr.f.eb,
                          rr.m.pe = rr.m.pe, rr.f.pe = rr.f.pe,
                          rr.m.pb = rr.m.pb, rr.f.pb = rr.f.pb,
                          rr.mf.bef = rr.mf.bef, rr.mf.exc = rr.mf.exc,
                          rr.mf.bef.cont = rr.mf.bef.cont, rr.mf.exc.cont = rr.mf.exc.cont)
        return(rrs)
    })}

## calculate the proportion of times a group of states occurred given a state.group.
prop.trans.arr <- function(seros, state.log, states, state.group) {
    apply(seros[state.log,,], 3, function(x) sum(x[,states]/x[,state.group])) / sum(state.log)
}
    

## Take seros as input and produce route of transmission breakdowns.
prop.trans <- function(seros, dat, browse = F) {
    ## Naming convention: pi = proportion of transmission, N(egative), P(positive), U(known) partner serostatus, A(live
    ## conditionality), b(efore couple duration), e(xtra-couple), p(artner transmitted),
    hh.log <- dat$ser==1
    mm.log <- dat$ser==2
    ff.log <- dat$ser==3
    ss.log <- dat$ser==4
    nhh <- sum(hh.log)
    nmm <- sum(mm.log)
    nff <- sum(ff.log)
    nss <- sum(ss.log)
    if(browse)  browser()
    ## male breakdown amongst observed M+F- couples, partner *N*egative
    pibNA  <- prop.trans.arr(seros, mm.log, c('mb.a1A', 'mb.a2A', 'mb.A'), 'mmA')
    pieNA  <- prop.trans.arr(seros, mm.log, c('me.a1A', 'me.a2A', 'me.A'), 'mmA')    
    ## female breakdown amongst observed M-F+ couples, partner *N*egative
    piNbA  <- prop.trans.arr(seros, ff.log, c('f.ba1A', 'f.ba2A', 'f.bA'), 'ffA')
    piNeA  <- prop.trans.arr(seros, ff.log, c('f.ea1A', 'f.ea2A', 'f.eA'), 'ffA')    
    ## male breakdown amongst observed M+F+ couples,  partner *P*ositive
    pibPA  <- prop.trans.arr(seros, hh.log, c('hb1b2A', 'hb2b1A', 'hbeA', 'hbpaA', 'hbpA'), 'hhA')
    piePA  <- prop.trans.arr(seros, hh.log, c('hebA', 'hepaA', 'hepA', 'he1e2A', 'he2e1A'), 'hhA')
    pipPA  <- prop.trans.arr(seros, hh.log, c('hpbaA', 'hpeaA', 'hpbA', 'hpeA'), 'hhA')
    pipPaA <- prop.trans.arr(seros, hh.log, c('hpbaA', 'hpeaA'), 'hhA') #acute
    pipPcA <- prop.trans.arr(seros, hh.log, c('hpbA', 'hpeA'),'hhA') #chronic
    ## female breakdown amongst observed M+F+ couples,  partner *P*ositive
    piPbA  <- prop.trans.arr(seros, hh.log, c('hb1b2A', 'hb2b1A', 'hebA', 'hpbaA', 'hpbA'), 'hhA')
    piPeA  <- prop.trans.arr(seros, hh.log, c('hbeA', 'hpeaA', 'hpeA', 'he1e2A', 'he2e1A'), 'hhA')
    piPpA  <- prop.trans.arr(seros, hh.log, c('hbpaA', 'hepaA', 'hbpA', 'hepA'), 'hhA')
    piPpaA <- prop.trans.arr(seros, hh.log, c('hbpaA', 'hepaA'), 'hhA') #acute
    piPpcA <- prop.trans.arr(seros, hh.log, c('hbpA', 'hepA'), 'hhA') #chronic
    ## male breakdown amongst infected males in any observed couples,  partner *U*nknown (bc could be either)
    pibUA  <- (pibNA*nmm + pibPA*nhh) / (nmm + nhh)
    pieUA  <- (pieNA*nmm + piePA*nhh) / (nmm + nhh)
    pipUA  <- (pipPA*nhh) / (nmm + nhh)
    pipUaA <- (pipPaA*nhh) / (nmm + nhh) #acute, keep in mind this is only proportion of all M+ acutely infected by their *current* partner
    pipUcA <- (pipPcA*nhh) / (nmm + nhh) #chronic, ditto above
    ## female breakdown amongst infected females in any observed couples,  partner *U*nknown (bc could be either)
    piUbA  <- (piNbA*nff + piPbA*nhh) / (nff + nhh)
    piUeA  <- (piNeA*nff + piPeA*nhh) / (nff + nhh)
    piUpA  <- (piPpA*nhh) / (nff + nhh)
    piUpaA <- (piPpaA*nhh) / (nff + nhh) #acute, see notes aboe
    piUpcA <- (piPpcA*nhh) / (nff + nhh) #chronic
    ## ####################################################################
    ## Route of transmission breakdowns for inferred
    ## pseudopopulation (unconditional on survival)
    ## ####################################################################
    ## Index infections, do with estimators summing over all
    ## infected couples and over all couples
    ## version 1 - all infected couples
    ## all infected and alive
    sro.names <- colnames(seros)
    all.infA.nms <- sro.names[(grepl("m", sro.names) | grepl("f", sro.names) | grepl("h", sro.names)) & grepl("A", sro.names)]
    all.infA.nms <- all.infA.nms[!all.infA.nms %in% c('mmA','ffA','hhA')]
    for(ii in 1:2) { ## with only infecteds & with everyone including S(usceptibles)
        if(ii==1)       cpls.to.use <- !ss.log
        if(ii==2)       {
            cpls.to.use <- rep(T, nrow(seros)) ## use all couples
            all.infA.nms <- c('s..', all.infA.nms) ## summing over alive states
        }
        mb1. <- apply(seros[cpls.to.use, c("mb.a1","mb.a2","mb.","hb1b2","hbe","hbpa","hbp"),], 3, rowSums)
        me1. <- apply(seros[cpls.to.use, c("me.a1","me.a2","me.","he1e2","hepa","hep"),], 3, rowSums)    
        ## female
        f.b1 <- apply(seros[cpls.to.use, c("f.ba1","f.ba2","f.b","hb2b1","heb","hpba","hpb"),], 3, rowSums)
        f.e1 <- apply(seros[cpls.to.use, c("f.ea1","f.ea2","f.e","he2e1","hpea","hpe"),], 3, rowSums)
        all.infA <- apply(seros[cpls.to.use, all.infA.nms,], 3, rowSums)
        ## number of inflated male before-couple index infections (in any couple)
        mb1.Infl <- colSums( mb1. / all.infA)
        ## number of inflated male extra-couple index infections (in any couple)
        me1.Infl <- colSums( me1. / all.infA)
        ## nufber of inflated male before-couple index infections (in any couple)
        f.b1Infl <- colSums( f.b1 / all.infA)
        ## number of inflated male extra-couple index infections (in any couple)
        f.e1Infl <- colSums( f.e1 / all.infA)
        ## number of inflated index infections (1 in each couple with an infection)
        IndInfl <- mb1.Infl + me1.Infl + f.b1Infl + f.e1Infl # note this includes h-classes..
        ## ####################################################################
        ## Proportion of index infections pooling across gender
        assign(paste0('piGb1.sumI',c('','S')[ii]), mb1.Infl / IndInfl)
        assign(paste0('piGe1.sumI',c('','S')[ii]), me1.Infl / IndInfl)
        assign(paste0('piG.b1sumI',c('','S')[ii]), f.b1Infl / IndInfl)
        assign(paste0('piG.e1sumI',c('','S')[ii]), f.e1Infl / IndInfl)
    }
    ## put them all in a dataframe
    pop.avs <- data.frame(
        ## conditional on survival
        pibNA, pieNA, # b/e in +- given A
        piNbA, piNeA, # b/e in -+ given A
        pibPA, piePA, pipPA, pipPaA, pipPcA, # b/e/p in male in ++ given A
        piPbA, piPeA, piPpA, piPpaA, piPpcA, # b/e/p in female in ++ given A
        pibUA, pieUA, pipUA, pipUaA, pipUcA, # b/e/p in male in any given A
        piUbA, piUeA, piUpA, piUpaA, piUpcA, # b/e/p in female in any given A
        ## unconditional on survival version 1
        piGb1.sumI, piGe1.sumI, # b/e index in males amongst all infected 
        piG.b1sumI, piG.e1sumI, # b/e index in females amongst all infected
        piGb1.sumIS, piGe1.sumIS, # b/e index in males amongst all infected 
        piG.b1sumIS, piG.e1sumIS) # b/e index in females amongst all infected
    return(pop.avs)
}



