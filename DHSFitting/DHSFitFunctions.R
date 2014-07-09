
pre.couple <- function(seros.active, pmb, pfb, pmb.a, pfb.a) {
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
    ## if(!fast) { ## if calculating probabilities of transmission with or without death (not used for fitting)
    ##     seros.active[,'mb.a1'] <- tp[,'s..'] * pmb * (1-pfb)
    ##     seros.active[,'mb.a2'] <- tp[,'mb.a1'] * (1 - pfb)
    ##     seros.active[,'mb.'] <- tp[,'mb.a2'] * (1 - pfb) + tp[,'mb.'] * (1 - pfb)
    ##     seros.active[,'f.ba1'] <- tp[,'s..'] * pfb * (1-pmb)
    ##     seros.active[,'f.ba2'] <- tp[,'f.ba1'] * (1 - pmb)
    ##     seros.active[,'f.b'] <- tp[,'f.ba2'] * (1 - pmb) + tp[,'f.b'] * (1 - pmb)
    ##     p.mfirst <- pmb / (pmb+pfb)
    ##     p.ffirst <- 1-p.mfirst
    ##     p.mfirst[is.na(p.mfirst)] <- 0
    ##     p.ffirst[is.na(p.ffirst)] <- 0                
    ##     seros.active[,'hb1b2'] <- tp[,'hb1b2'] + p.mfirst  *  tp[,'s..'] * pmb * pfb + (tp[,'mb.a1'] + tp[,'mb.a2'] + tp[,'mb.'])  *  pfb
    ##     seros.active[,'hb2b1'] <- tp[,'hb2b1'] + p.ffirst  *  tp[,'s..'] * pmb * pfb + (tp[,'f.ba1'] + tp[,'f.ba2'] + tp[,'f.b'])  *  pmb
    ## }
    if(sum(is.na(seros.active))>0) browser()
    return(seros.active)
}

within.couple <- function(seros.active, pme, pfe, pmp, pfp, pmp.ac, pfp.ac, pme.a, pfe.a, pmp.a, pfp.a, pmp.a.ac, pfp.a.ac) {
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
    ## if(!fast) { ## if calculating probabilities of transmission with or without death (not used for fitting)
    ##     seros.active[,'mb.a1'] <- 0 
    ##     seros.active[,'mb.a2'] <- tp[,'mb.a1'] * (1-pfe)*(1-pfp.ac)
    ##     seros.active[,'mb.']   <- tp[,'mb.a2'] * (1-pfe)*(1-pfp.ac) + tp[,'mb.'] * (1-pfe)*(1-pfp)            
    ##     seros.active[,'me.a1'] <- tp[,'s..'] * pme*(1-pfe)
    ##     seros.active[,'me.a2'] <- tp[,'me.a1']*(1-pfe)*(1-pfp.ac) 
    ##     seros.active[,'me.']   <- tp[,'me.a2']*(1-pfe)*(1-pfp.ac) + tp[,'me.']*(1-pfe)*(1-pfp) 
    ##     seros.active[,'f.ba1'] <- 0
    ##     seros.active[,'f.ba2'] <- tp[,'f.ba1']*(1-pme)*(1-pmp.ac)
    ##     seros.active[,'f.b']   <- tp[,'f.ba2']*(1-pme)*(1-pmp.ac) + tp[,'f.b']*(1-pme)*(1-pmp)            
    ##     seros.active[,'f.ea1'] <- tp[,'s..']*pfe*(1-pme)
    ##     seros.active[,'f.ea2'] <- tp[,'f.ea1']*(1-pme)*(1-pmp.ac)
    ##     seros.active[,'f.e']   <- tp[,'f.ea2']*(1-pme)*(1-pmp.ac) + tp[,'f.e']*(1-pme)*(1-pmp)
    ##     ##          hb1b2'] <- hb1b2'] # Doesn't change during couple duration
    ##     ##          hb2b1'] <- hb2b1'] # Doesn't change during couple duration                
    ##     seros.active[,'hbe']  <- tp[,'hbe'] + (tp[,'mb.a1']+ tp[,'mb.a2'])*(1-pfp.ac)*pfe + tp[,'mb.']*(1-pfp)*pfe
    ##     seros.active[,'heb']  <- tp[,'heb'] + (tp[,'f.ba1']+ tp[,'f.ba2'])*(1-pmp.ac)*pme + tp[,'f.b']*(1-pmp)*pme
    ##     seros.active[,'hbpa'] <- tp[,'hbpa']+ (tp[,'mb.a1']+ tp[,'mb.a2'])*pfp.ac
    ##     seros.active[,'hbp']  <- tp[,'hbp'] + tp[,'mb.']*pfp            
    ##     seros.active[,'hpba'] <- tp[,'hpba']+ (tp[,'f.ba1']+ tp[,'f.ba2'])*pmp.ac
    ##     seros.active[,'hpb']  <- tp[,'hpb'] + tp[,'f.b']*pmp
    ##     seros.active[,'hepa'] <- tp[,'hepa']+ (tp[,'me.a1']+ tp[,'me.a2'])*pfp.ac
    ##     seros.active[,'hep']  <- tp[,'hep'] + tp[,'me.']*pfp            
    ##     seros.active[,'hpea'] <- tp[,'hpea']+ (tp[,'f.ea1']+ tp[,'f.ea2'])*pmp.ac
    ##     seros.active[,'hpe']  <- tp[,'hpe'] + tp[,'f.e']*pmp            
    ##     ## for individuals infected in the same month, assign
    ##     ## the order of infection based on competing risks
    ##     ## formula, but if the denominator is 0, replace both
    ##     ## with 0 to avoid errors.
    ##     p.mfirst <- pme / (pme+pfe)
    ##     p.ffirst <- 1-p.mfirst
    ##     p.mfirst[is.na(p.mfirst)] <- 0
    ##     p.ffirst[is.na(p.ffirst)] <- 0                
    ##     seros.active[,'he1e2'] <- tp[,'he1e2'] + p.mfirst * tp[,'s..']*pme*pfe + (tp[,'me.a1'] + tp[,'me.a2'])*(1-pfp.ac)*pfe + tp[,'me.']*(1-pfp)*pfe
    ##     seros.active[,'he2e1'] <- tp[,'he2e1'] + p.ffirst * tp[,'s..']*pme*pfe + (tp[,'f.ea1'] + tp[,'f.ea2'])*(1-pmp.ac)*pme + tp[,'f.e']*(1-pmp)*pme
    ## }
    return(seros.active)
}

## get logical matrices for whether male, female, or either partner is sexually active for each
## month and country-level prevalence for each couple's tt-th month of before-couple duration. For
## use in pre-couple transmission function.
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
## get logical matrices for whether the couple is formed yet and country-level prevalence for each
## couple's tt-th month of couple duration. Also survival & ART coverage. For use in within-couple
## transmission function.
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

pcalc <- function(pars, dat, browse = F, 
                  give.ser = F, ## return serostatus probability matrix
                  acute.sc = 7, ## acute phase infectiousness
                  ## using global variables for these
                  ## partner.arv = F, ## ART coverage affects within-partnership transmission?
                  ## low.coverage.arv = F, ## if T, only 50% of those on ART are not infectious
                  lrho.sd = 1/2) ## sd on lrho prior
    {
        if(browse) browser()
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
            state.var.nms <- c('s..',
                               'mb.a1', 'mb.a2', 'mb.', 'me.a1', 'me.a2', 'me.', #M+F- by routes & acute phase (a)
                               'f.ba1', 'f.ba2', 'f.b', 'f.ea1', 'f.ea2', 'f.e', #M-F+ by routes & acute phase (a)
                               'hb1b2', 'hb2b1', 'hbe', 'heb', ## rest are all M+F+
                               'hbpa', 'hpba', 'hepa', 'hpea', ## acute infected by partner
                               'hbp', 'hpb', 'hep', 'hpe', ## chronic infected by partner
                               'he2e1', 'he1e2') ## both extra-couply infected, different orders
            state.var.nms <- c(state.var.nms, paste0(state.var.nms[-1],'A')) ## joint with alive at the end
            ## #####################################################################################          
            ## i.e., hbpa is a ++ couple in which the male was inf *b*efore
            ## couple formation & the female by her *p*artner while he was *a*cutely infectious
            ## #####################################################################################
            seros <- matrix(0, K, length(state.var.nms), dimnames = list(NULL,state.var.nms))
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
                seros[active,] <- pre.couple(seros[active,], pmb=pmb, pfb=pfb, pmb.a=pmb.a, pfb.a=pfb.a) ## update serostates
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
                s.p.f <- within.msurv[ff,tt]
                pmp.a <- pmp * s.p.m                ## Transmission probabilities from partner (jointly with survival)
                pfp.a <- pfp * s.p.f
                pmp.a.ac <- pmp.ac * s.p.m                ## acute from partner (jointly with survival)
                pfp.a.ac <- pfp.ac * s.p.f
                pme.a <- pme * s.p.m                ## extra-couple
                pfe.a <- pfe * s.p.f
                ## Update serostates for active couples
                seros[ff,] <- within.couple(seros[ff,], pme=pme, pfe=pfe, pmp=pmp, pfp=pfp, pmp.ac=pmp.ac, pfp.ac = pfp.ac,
                                            pme.a = pme.a, pfe.a = pfe.a, pmp.a = pmp.a, pfp.a = pfp.a, pmp.a.ac = pmp.a.ac, pfp.a.ac = pfp.a.ac)
            }
            ## sum.nms <- c('ss','mm','ff','hh','mmA','ffA','hhA') ## summary serostatuses (observable) unconditional on survival
            sum.nms <- c('ss','mmA','ffA','hhA') ## summary serostatuses (observable)
            sero.sums <- matrix(NA, K, length(sum.nms), dimnames = list(NULL,sum.nms))
            sero.sums[,'ss'] <- seros[,'s..']
            ## sero.sums[,'mm'] <- rowSums(seros[,c('mb.a1', 'me.a1', 'mb.a2', 'me.a2', 'mb.', 'me.')])
            ## sero.sums[,'ff'] <- rowSums(seros[,c('f.ba1', 'f.ea1', 'f.ba2', 'f.ea2', 'f.b', 'f.e')])
            ## sero.sums[,'hh'] <- rowSums(seros[,c('hb1b2', 'hb2b1', 'hbe', 'heb', 'hbpa', 'hpba', 'hepa', 'hpea', 'hbp', 'hpb', 'hep', 'hpe', 'he1e2', 'he2e1')])
            sero.sums[,'mmA'] <- rowSums(seros[,c('mb.a1A', 'me.a1A', 'mb.a2A', 'me.a2A', 'mb.A', 'me.A')])
            sero.sums[,'ffA'] <- rowSums(seros[,c('f.ba1A', 'f.ea1A', 'f.ba2A', 'f.ea2A', 'f.bA', 'f.eA')])
            sero.sums[,'hhA'] <- rowSums(seros[,c('hb1b2A', 'hb2b1A', 'hbeA', 'hebA', 'hbpaA', 'hpbaA', 'hepaA', 'hpeaA', 'hbpA', 'hpbA', 'hepA', 'hpeA', 'he1e2A', 'he2e1A')])
            ## if(survive) { ## must NORMALIZE probabilities to 1 for likelihood!
            pser.a <- sero.sums[,c('hhA', 'mmA', 'ffA', 'ss')]
            probs <- pser.a[cbind(1:K,dat$ser)] / rowSums(pser.a) # accounting for survival
            ## }else{
            ##     pser <- sero.sums[,c('hh', 'mm', 'ff', 'ss')]
            ##     probs <- pser[cbind(1:K,dat$ser)] ## if not accounting for survival
            ## }
            if(sum(probs==0)==0) { # if non of the serotatuses occur with 0 probability in the current model
                lprob <- sum(log(probs)) + dnorm(log(rho), log(trans.ratio), lrho.sd, log = T) ## add llikelihood and lprior
            }else{ # if some of the serostatuses are 0, then the current parameters have 0 probability
                lprob <- -Inf
            }
        }
        return(lprob)
    }


##         if(trace & sim) # calculate expected route of transmission breakdowns for couples with unknown (or simulated serostatus)
##           {
## ######################################################################
##             ## Route of transmission breakdowns for observed couples
##             ## (conditional on survival)
## ######################################################################
##             ## male breakdown amongst observed M+F- couples, partner *N*egative
##             pibNA <- mean((mb.a1A + mb.a2A + mb.A) / mmA, na.rm = T) 
##             pieNA <- mean((me.a1A + me.a2A + me.A) / mmA, na.rm = T) 
##             ## female breakdown amongst observed M-F+ couples, partner *N*egative
##             piNbA <- mean((f.ba1A + f.ba2A + f.bA) / ffA, na.rm = T)
##             piNeA <- mean((f.ea1A + f.ea2A + f.eA) / ffA, na.rm = T)
##             ## male breakdown amongst observed M+F+ couples,  partner *P*ositive
##             pibPA <- mean((hb1b2A + hb2b1A + hbeA + hbpaA + hbpA) / hhA, na.rm = T) 
##             piePA <- mean((hebA + hepaA + hepA + he1e2A + he2e1A) / hhA, na.rm = T) 
##             pipPA <- mean((hpbaA + hpeaA + hpbA + hpeA) / hhA, na.rm = T)
##             pipPaA <- mean((hpbaA + hpeaA) / hhA, na.rm = T) # males infected by their partner in M+F+ couples during her acute phase.
##             pipPcA <- mean((hpbA + hpeA) / hhA, na.rm = T) # chronic phase
##             ## female breakdown amongst observed M+F+ couples,  partner *P*ositive
##             piPbA <- mean((hb1b2A + hb2b1A + hebA + hpbaA + hpbA) / hhA, na.rm = T) 
##             piPeA <- mean((hbeA + hpeaA + hpeA + he1e2A + he2e1A) / hhA, na.rm = T) 
##             piPpA <- mean((hbpaA + hepaA + hbpA + hepA) / hhA, na.rm = T)
##             piPpaA <- mean((hbpaA + hepaA) / hhA, na.rm = T) # acute phase
##             piPpcA <- mean((hbpA + hepA) / hhA, na.rm = T) # chronic
##             ## male breakdown amongst infected males in any observed couples,  partner *U*nknown (bc could be either)
##             pibUA <- mean((mb.a1A + mb.a2A + mb.A + hb1b2A + hb2b1A + hbeA + hbpaA + hbpA) / (mmA + hhA), na.rm = T)
##             pieUA <- mean((me.a1A + me.a2A + me.A + hebA + hepaA + hepA + he1e2A + he2e1A) / (mmA + hhA), na.rm = T)
##             pipUA <- mean((hpbaA + hpeaA + hpbA + hpeA) / (mmA + hhA), na.rm = T)
##             pipUaA <- mean((hpbaA + hpeaA) / (mmA + hhA), na.rm = T) #acute
##             pipUcA <- mean((hpbA + hpeA) / (mmA + hhA), na.rm = T) #chronic
##             ## female breakdown amongst infected females in any observed couples,  partner *U*nknown (bc could be either)
##             piUbA <- mean((f.ba1A + f.ba2A + f.bA + hb1b2A + hb2b1A + hebA + hpbaA + hpbA) / (ffA + hhA), na.rm = T)
##             piUeA <- mean((f.ea1A + f.ea2A + f.eA + hbeA + hpeaA + hpeA + he1e2A + he2e1A) / (ffA + hhA), na.rm = T)
##             piUpA <- mean((hbpaA + hepaA + hbpA + hepA) / (ffA + hhA), na.rm = T)
##             piUpaA <- mean((hbpaA + hepaA) / (ffA + hhA), na.rm = T) #acute
##             piUpcA <- mean((hbpA + hepA) / (ffA + hhA), na.rm = T) #chronic
## ######################################################################
##             pop.avs <- data.frame(
##                                   ## conditional on survival
##                                   pibNA, pieNA, # b/e in +- given A
##                                   piNbA, piNeA, # b/e in -+ given A
##                                   pibPA, piePA, pipPA, pipPaA, pipPcA, # b/e/p in male in ++ given A
##                                   piPbA, piPeA, piPpA, piPpaA, piPpcA, # b/e/p in female in ++ given A
##                                   pibUA, pieUA, pipUA, pipUaA, pipUcA, # b/e/p in male in any given A
##                                   piUbA, piUeA, piUpA, piUpaA, piUpcA) # b/e/p in female in any given A
##           }
##         if(trace & !sim)
##           {
##             ######################################################################
##             ## Route of transmission breakdowns for observed couples
##             ## (conditional on survival)
##             ######################################################################
##             ## male breakdown amongst observed M+F- couples, partner *N*egative
##             pibNA <- sum((mb.a1A[mm.log] + mb.a2A[mm.log] + mb.A[mm.log]) / mmA[mm.log]) / sum(mm.log)
##             pieNA <- sum((me.a1A[mm.log] + me.a2A[mm.log] + me.A[mm.log]) / mmA[mm.log]) / sum(mm.log)
##             ## female breakdown amongst observed M-F+ couples, partner *N*egative
##             piNbA <- sum((f.ba1A[ff.log] + f.ba2A[ff.log] + f.bA[ff.log]) / ffA[ff.log]) / sum(ff.log)
##             piNeA <- sum((f.ea1A[ff.log] + f.ea2A[ff.log] + f.eA[ff.log]) / ffA[ff.log]) / sum(ff.log)
##             ## male breakdown amongst observed M+F+ couples,  partner *P*ositive
##             pibPA <- sum((hb1b2A[hh.log] + hb2b1A[hh.log] + hbeA[hh.log] + hbpaA[hh.log] + hbpA[hh.log]) / hhA[hh.log]) / sum(hh.log)
##             piePA <- sum((hebA[hh.log] + hepaA[hh.log] + hepA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]) / sum(hh.log)
##             pipPA <- sum((hpbaA[hh.log] + hpeaA[hh.log] + hpbA[hh.log] + hpeA[hh.log]) / hhA[hh.log]) / sum(hh.log)
##             pipPaA <- sum((hpbaA[hh.log] + hpeaA[hh.log]) / hhA[hh.log]) / sum(hh.log) #acute
##             pipPcA <- sum((hpbA[hh.log] + hpeA[hh.log]) / hhA[hh.log]) / sum(hh.log) #chronic
##             ## female breakdown amongst observed M+F+ couples,  partner *P*ositive
##             piPbA <- sum((hb1b2A[hh.log] + hb2b1A[hh.log] + hebA[hh.log] + hpbaA[hh.log] + hpbA[hh.log]) / hhA[hh.log]) / sum(hh.log)
##             piPeA <- sum((hbeA[hh.log] + hpeaA[hh.log] + hpeA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]) / sum(hh.log)
##             piPpA <- sum((hbpaA[hh.log] + hepaA[hh.log] + hbpA[hh.log] + hepA[hh.log]) / hhA[hh.log]) / sum(hh.log)
##             piPpaA <- sum((hbpaA[hh.log] + hepaA[hh.log]) / hhA[hh.log]) / sum(hh.log) #acute
##             piPpcA <- sum((hbpA[hh.log] + hepA[hh.log]) / hhA[hh.log]) / sum(hh.log) #chronic
##             ## male breakdown amongst infected males in any observed couples,  partner *U*nknown (bc could be either)
##             pibUA <- (pibNA*sum(mm.log) + pibPA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
##             pieUA <- (pieNA*sum(mm.log) + piePA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
##             pipUA <- (pipPA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
##             pipUaA <- (pipPaA*sum(hh.log)) / (sum(mm.log) + sum(hh.log)) #acute, keep in mind this is only proportion of all M+ acutely infected by their *current* partner
##             pipUcA <- (pipPcA*sum(hh.log)) / (sum(mm.log) + sum(hh.log)) #chronic, ditto above
##             ## female breakdown amongst infected females in any observed couples,  partner *U*nknown (bc could be either)
##             piUbA <- (piNbA*sum(ff.log) + piPbA*sum(hh.log)) / (sum(ff.log) + sum(hh.log))
##             piUeA <- (piNeA*sum(ff.log) + piPeA*sum(hh.log)) / (sum(ff.log) + sum(hh.log))
##             piUpA <- (piPpA*sum(hh.log)) / (sum(ff.log) + sum(hh.log))
##             piUpaA <- (piPpaA*sum(hh.log)) / (sum(ff.log) + sum(hh.log)) #acute, see notes aboe
##             piUpcA <- (piPpcA*sum(hh.log)) / (sum(ff.log) + sum(hh.log)) #chronic
##             ######################################################################
##             ## Give pieNA, piNeA, piePA, piPeA, for each couple
##             if(give.pis)
##               {
##                 ## probability infection was extracouple given ser
##                 piCe.A <- rep(NA, K)
##                 piC.eA <- rep(NA, K)
##                 piCe.A[mm.log] <- (me.a1A[mm.log] + me.a2A[mm.log] + me.A[mm.log]) / mmA[mm.log]
##                 piCe.A[hh.log] <- (hebA[hh.log] + hepaA[hh.log] + hepA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]
##                 piC.eA[ff.log] <- (f.ea1A[ff.log] + f.ea2A[ff.log] + f.eA[ff.log]) / ffA[ff.log]
##                 piC.eA[hh.log] <- (hbeA[hh.log] + hpeaA[hh.log] + hpeA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]
##                 pis <- data.frame(piCe.A, piC.eA)
##               }
##             ######################################################################
##             ## Route of transmission breakdowns for inferred
##             ## pseudopopulation (unconditional on survival)
##             ######################################################################
##             ######################################################################
##             ## Index infections, do with estimators summing over all
##             ## infected couples and over all couples
##             ## version 1 - all infected couples
##             mb1. <- rowSums(allstates[!ss.log, c("mb.a1","mb.a2","mb.","hb1b2","hbe","hbpa","hbp")])
##             me1. <- rowSums(allstates[!ss.log, c("me.a1","me.a2","me.","he1e2","hepa","hep")])
##             ## female
##             f.b1 <- rowSums(allstates[!ss.log, c("f.ba1","f.ba2","f.b","hb2b1","heb","hpba","hpb")])
##             f.e1 <- rowSums(allstates[!ss.log, c("f.ea1","f.ea2","f.e","he2e1","hpea","hpe")])
##             all.infA <- rowSums(allstates[!ss.log,names(allstates)[(grepl("m", names(allstates)) | grepl("f", names(allstates)) | grepl("h", names(allstates))) & grepl("A", names(allstates))]])
##             ## number of inflated male before-couple index infections (in any couple)
##             mb1.Infl <- sum( mb1. / all.infA)
##             ## number of inflated male extra-couple index infections (in any couple)
##             me1.Infl <- sum( me1. / all.infA)
##             ## nufber of inflated male before-couple index infections (in any couple)
##             f.b1Infl <- sum( f.b1 / all.infA)
##             ## number of inflated male extra-couple index infections (in any couple)
##             f.e1Infl <- sum( f.e1 / all.infA)
##             ## number of inflated index infections (1 in each couple with an infection)
##             IndInfl <- mb1.Infl + me1.Infl + f.b1Infl + f.e1Infl # note this includes h-classes..
##             ######################################################################
##             ## Proportion of index infections pooling across gender
##             piGb1.sumI <- mb1.Infl / IndInfl
##             piGe1.sumI <- me1.Infl / IndInfl
##             piG.b1sumI <- f.b1Infl / IndInfl
##             piG.e1sumI <- f.e1Infl / IndInfl
##             ## ## version 2 - sum over all couples
##             mb1. <- rowSums(allstates[, c("mb.a1","mb.a2","mb.","hb1b2","hbe","hbpa","hbp")])
##             me1. <- rowSums(allstates[, c("me.a1","me.a2","me.","he1e2","hepa","hep")])
##             ## female
##             f.b1 <- rowSums(allstates[, c("f.ba1","f.ba2","f.b","hb2b1","heb","hpba","hpb")])
##             f.e1 <- rowSums(allstates[, c("f.ea1","f.ea2","f.e","he2e1","hpea","hpe")])
##             all.infA <- rowSums(allstates[,c("s..", names(allstates)[(grepl("m", names(allstates)) | grepl("f", names(allstates)) | grepl("h", names(allstates))) & grepl("A", names(allstates))])])
##             ## number of inflated male before-couple index infections (in any couple)
##             mb1.Infl <- sum( mb1. / all.infA)
##             ## number of inflated male extra-couple index infections (in any couple)
##             me1.Infl <- sum( me1. / all.infA)
##             ## nufber of inflated male before-couple index infections (in any couple)
##             f.b1Infl <- sum( f.b1 / all.infA)
##             ## number of inflated male extra-couple index infections (in any couple)
##             f.e1Infl <- sum( f.e1 / all.infA)
##             ## number of inflated index infections (1 in each couple with an infection)
##             IndInfl <- mb1.Infl + me1.Infl + f.b1Infl + f.e1Infl
## ######################################################################
##             ## Proportion of index infections pooling across gender
##             piGb1.sumIS <- mb1.Infl / IndInfl
##             piGe1.sumIS <- me1.Infl / IndInfl
##             piG.b1sumIS <- f.b1Infl / IndInfl
##             piG.e1sumIS <- f.e1Infl / IndInfl
##             ## put them all in a dataframe
##             pop.avs <- data.frame(
##                                   ## conditional on survival
##                                   pibNA, pieNA, # b/e in +- given A
##                                   piNbA, piNeA, # b/e in -+ given A
##                                   pibPA, piePA, pipPA, pipPaA, pipPcA, # b/e/p in male in ++ given A
##                                   piPbA, piPeA, piPpA, piPpaA, piPpcA, # b/e/p in female in ++ given A
##                                   pibUA, pieUA, pipUA, pipUaA, pipUcA, # b/e/p in male in any given A
##                                   piUbA, piUeA, piUpA, piUpaA, piUpcA, # b/e/p in female in any given A
##                                   ## unconditional on survival version 1
##                                   piGb1.sumI, piGe1.sumI, # b/e index in males amongst all infected 
##                                   piG.b1sumI, piG.e1sumI, # b/e index in females amongst all infected
##                                   piGb1.sumIS, piGe1.sumIS, # b/e index in males amongst all infected 
##                                   piG.b1sumIS, piG.e1sumIS) # b/e index in females amongst all infected
##             ######################################################################
##             ## Project incidence forward 12 months for each of the 3
##             ## couple types (ss, mm, ff) for each country in the data
##             ## set (because they have different population prevalences
##             num.country <- length(unique(dat$epic.ind))
##             cc.inds <- unique(dat$epic.ind)
##             ## concordant negative
##             ss12.ssL <- rep(1, num.country)
##             mm12a1.ssL <- rep(0, num.country)
##             mm12a2.ssL <- rep(0, num.country)
##             mm12.ssL <- rep(0, num.country)            
##             ff12a1.ssL <- rep(0, num.country)
##             ff12a2.ssL <- rep(0, num.country)
##             ff12.ssL <- rep(0, num.country)            
##             hh12.ssL <- rep(0, num.country)            
##             ## male positive discordant
##             mm12a1.mmL <- rep(1, num.country)
##             mm12a2.mmL <- rep(1, num.country)
##             mm12.mmL <- rep(1, num.country)            
##             hh12.mmL <- rep(0, num.country)            
##             ## female positive discordant
##             ff12a1.ffL <- rep(1, num.country)
##             ff12a2.ffL <- rep(1, num.country)
##             ff12.ffL <- rep(1, num.country)            
##             hh12.ffL <- rep(0, num.country)
##             ## initialize pis
##             pi.m.part12.ss.ac <- 0      #acute
##             pi.f.part12.ss.ac <- 0      #acute
##             pi.m.part12.ss <- 0
##             pi.f.part12.ss <- 0
##             pi.m.exc12.ss <- 0
##             pi.f.exc12.ss <- 0
##             pi.f.part12.mm.ac <- 0      #acute
##             pi.f.part12.mm <- 0            
##             pi.f.exc12.mm <- 0            
##             pi.m.part12.ff.ac <- 0      #acute
##             pi.m.part12.ff <- 0            
##             pi.m.exc12.ff <- 0            
##             for(tt in 1:12)
##               {
##                 if(partner.arv)         # put partner's on ART with some probability assigned related to ART coverage
##                   {
##                     p.m.part <- 1 - exp(-bmp * (1 - cov.scalar * art.prev[1332+tt-1, cc.inds]))
##                     p.f.part <- 1 - exp(-bfp * (1 - cov.scalar * art.prev[1332+tt-1, cc.inds]))
##                     p.m.part.ac <- 1 - exp(-acute.sc * bmp * (1 - cov.scalar * art.prev[1332+tt-1, cc.inds])) # acute
##                     p.f.part.ac <- 1 - exp(-acute.sc * bfp * (1 - cov.scalar * art.prev[1332+tt-1, cc.inds])) # acute
##                   }
##                 ######################################################################
##                 ## Transmission probabilities
##                 ## probability infected extracouply in various months of 2011
##                 p.m.exc <- 1 - exp(-bme*epicf[1332+tt-1, cc.inds])
##                 p.f.exc <- 1 - exp(-bfe*epicm[1332+tt-1, cc.inds])
##                 ## concordant negative couples
##                 ss12.ss   <- ss12.ssL*(1-p.m.exc)*(1-p.f.exc)
##                 mm12a1.ss <- ss12.ssL*p.m.exc*(1-p.f.exc)
##                 mm12a2.ss <- mm12a1.ssL*(1-p.f.exc)*(1-p.f.part.ac)
##                 mm12.ss   <- mm12a2.ssL*(1-p.f.exc)*(1-p.f.part.ac) + mm12.ssL*(1-p.f.exc)*(1-p.f.part)
##                 ff12a1.ss <- ss12.ssL*p.f.exc*(1-p.m.exc)
##                 ff12a2.ss <- ff12a1.ssL*(1-p.m.exc)*(1-p.m.part.ac)
##                 ff12.ss   <- ff12a2.ssL*(1-p.m.exc)*(1-p.m.part.ac) + ff12.ssL*(1-p.m.exc)*(1-p.m.part)
##                 hh12.ss <- hh12.ssL + ss12.ssL* p.m.exc*p.f.exc +
##                     (mm12a1.ssL + mm12a2.ssL)*(p.f.part.ac + (1-p.f.part.ac)*p.f.exc) + mm12.ssL*(p.f.part + (1-p.f.part)*p.f.exc) +
##                     (ff12a1.ssL + ff12a2.ssL)*(p.m.part.ac + (1-p.m.part.ac)*p.m.exc) + ff12.ssL*(p.m.part + (1-p.m.part)*p.m.exc)
##                 pi.m.part12.ss <- pi.m.part12.ss + (ff12a1.ssL + ff12a2.ssL)*p.m.part.ac + ff12.ssL*p.m.part
##                 pi.f.part12.ss <- pi.f.part12.ss + (mm12a1.ssL + mm12a2.ssL)*p.f.part.ac + mm12.ssL*p.f.part        
##                 pi.m.part12.ss.ac <- pi.m.part12.ss.ac + (ff12a1.ssL + ff12a2.ssL)*p.m.part.ac  # acute
##                 pi.f.part12.ss.ac <- pi.f.part12.ss.ac + (mm12a1.ssL + mm12a2.ssL)*p.f.part.ac  # acute
##                 pi.m.exc12.ss <- pi.m.exc12.ss + (ss12.ssL + (ff12a1.ssL + ff12a2.ssL)*(1-p.m.part.ac) + ff12.ssL*(1-p.m.part))*p.m.exc
##                 pi.f.exc12.ss <- pi.f.exc12.ss + (ss12.ssL + (mm12a1.ssL + mm12a2.ssL)*(1-p.f.part.ac) + mm12.ssL*(1-p.f.part))*p.f.exc
##                 ## male positive couples & female seroconversion,
##                 ##  assume there are no prevalent serodiscordant couples with individuals in acute stage at DHS ****
##                 mm12a1.mm <- 0          # all M+ couples  are already M+ at start of 12 month projection
##                 mm12a2.mm <- 0
##                 mm12.mm   <- mm12.mmL*(1-p.f.exc)*(1-p.f.part)                
##                 hh12.mm   <- hh12.mmL + mm12.mmL*(p.f.part + (1-p.f.part)*p.f.exc) 
##                 pi.f.part12.mm <- pi.f.part12.mm + mm12.mmL*p.f.part        
##                 pi.f.exc12.mm <- pi.f.exc12.mm + mm12.mmL*(1-p.f.part)*p.f.exc
##                 ## female positive couples & male seroconversion                  
##                 ff12a1.ff <- 0
##                 ff12a2.ff <- 0
##                 ff12.ff   <- ff12.ffL*(1-p.m.exc)*(1-p.m.part)                
##                 hh12.ff   <- hh12.ffL + ff12.ffL*(p.m.part + (1-p.m.part)*p.m.exc)
##                 pi.m.part12.ff <- pi.m.part12.ff + ff12.ffL*p.m.part
##                 pi.m.exc12.ff <- pi.m.exc12.ff + ff12.ffL*(1-p.m.part)*p.m.exc
##                 ss12.ssL   <- ss12.ss
##                 mm12a1.ssL <- mm12a1.ss
##                 mm12a2.ssL <- mm12a2.ss
##                 mm12.ssL   <- mm12.ss                
##                 ff12a1.ssL <- ff12a1.ss
##                 ff12a2.ssL <- ff12a2.ss
##                 ff12.ssL   <- ff12.ss                
##                 hh12.ssL   <- hh12.ss
##                 ## male positive discordant
##                 mm12a1.mmL <- mm12a1.mm
##                 mm12a2.mmL <- mm12a2.mm
##                 mm12.mmL   <- mm12.mm
##                 hh12.mmL   <- hh12.mm
##                 ## female positive discordant
##                 ff12a1.ffL <- ff12a1.ff
##                 ff12a2.ffL <- ff12a2.ff
##                 ff12.ffL   <- ff12.ff                
##                 hh12.ffL   <- hh12.ff
##               }
##             n.m.part.dc <- 0
##             n.m.part.cc.ac <- 0
##             n.m.part.cc <- 0            
##             n.f.part.dc <- 0
##             n.f.part.cc.ac <- 0
##             n.f.part.cc <- 0            
##             n.m.exc.dc <- 0
##             n.m.exc.cc <- 0
##             n.f.exc.dc <- 0
##             n.f.exc.cc <- 0
##             ## add all the different countries incidence by scaling by serotype
##             for(cc in 1:num.country)
##               {
##                 n.m.part.dc <- n.m.part.dc + pi.m.part12.ff[cc]*sum(ff.log & dat$epic.ind == cc.inds[cc])
##                 n.f.part.dc <- n.f.part.dc + pi.f.part12.mm[cc]*sum(mm.log & dat$epic.ind == cc.inds[cc])                
##                 n.m.part.cc <- n.m.part.cc + pi.m.part12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])
##                 n.f.part.cc <- n.f.part.cc + pi.f.part12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])                
##                 n.m.part.cc.ac <- n.m.part.cc.ac + pi.m.part12.ss.ac[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])
##                 n.f.part.cc.ac <- n.f.part.cc.ac + pi.f.part12.ss.ac[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])                
##                 n.m.exc.dc <- n.m.exc.dc + pi.m.exc12.ff[cc]*sum(ff.log & dat$epic.ind == cc.inds[cc])
##                 n.f.exc.dc <- n.f.exc.dc + pi.f.exc12.mm[cc]*sum(mm.log & dat$epic.ind == cc.inds[cc])                
##                 n.m.exc.cc <- n.m.exc.cc + pi.m.exc12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])
##                 n.f.exc.cc <- n.f.exc.cc + pi.f.exc12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])                
##               }
##             n.m.dc <- n.m.part.dc + n.m.exc.dc
##             n.f.dc <- n.f.part.dc + n.f.exc.dc
##             n.m.part.tot <- n.m.part.dc + n.m.part.cc
##             n.f.part.tot <- n.f.part.dc + n.f.part.cc            
##             n.m.exc.tot <- n.m.exc.dc + n.m.exc.cc
##             n.f.exc.tot <- n.f.exc.dc + n.f.exc.cc
##             proj12 <- data.frame(n.m.part.dc, n.f.part.dc, # incidence per 1000
##                                  n.m.part.cc, n.f.part.cc,
##                                  n.m.part.cc.ac, n.f.part.cc.ac,
##                                  n.m.exc.dc, n.f.exc.dc,
##                                  n.m.exc.cc, n.f.exc.cc,
##                                  n.m.part.tot, n.f.part.tot,
##                                  n.m.exc.tot, n.f.exc.tot) / sum(!hh.log) * 1000
##             prop.exc.m <- n.m.exc.tot / (n.m.exc.tot + n.m.part.tot)
##             prop.exc.f <- n.f.exc.tot / (n.f.exc.tot + n.f.part.tot)
##             prop.exc.m.dc <- n.m.exc.dc / n.m.dc
##             prop.exc.f.dc <- n.f.exc.dc / n.f.dc
##             proj12 <- data.frame(proj12, prop.exc.m, prop.exc.f, prop.exc.m.dc, prop.exc.f.dc) 
##             ## relative rate of transmission coefficient extracouply vs before relationship
##             rr.m.eb <- bme/bmb
##             rr.f.eb <- bfe/bfb
##             rr.m.pe <- bmp/bme
##             rr.f.pe <- bfp/bfe
##             rr.m.pb <- bmp/bmb        #partner to before
##             rr.f.pb <- bfp/bfb
##             ## relative rate of transmission coefficient extracouply and before relationship between males and females
##             rr.mf.bef <- bmb/bfb
##             rr.mf.exc <- bme/bfe
##             ## rho is the last one
##             ## relative rate of contact/risk paramter (i.e. accounting
##             ## for difference in per coital act probability as estimated
##             ## from within partnership transmission.
##             rr.mf.bef.cont <- rr.mf.bef * rho
##             rr.mf.exc.cont <- rr.mf.exc * rho
##             rrs <- data.frame(rr.m.eb = rr.m.eb, rr.f.eb = rr.f.eb,
##                               rr.m.pe = rr.m.pe, rr.f.pe = rr.f.pe,
##                               rr.m.pb = rr.m.pb, rr.f.pb = rr.f.pb,
##                               rr.mf.bef = rr.mf.bef, rr.mf.exc = rr.mf.exc,
##                               rr.mf.bef.cont = rr.mf.bef.cont, rr.mf.exc.cont = rr.mf.exc.cont)
##           }
##         if(sim) {# if simulating data
##             probs <- NA
##             lprob <- NA
##             ## create couple state probability *A*live & *D*ead
##             sim.probs <- data.frame(s..A = s..,
##                                     mb.a1A,   mb.a1D   =  mb.a1 - mb.a1A,
##                                     me.a1A,   me.a1D   =  me.a1 - me.a1A,
##                                     f.ba1A,   f.ba1D   =  f.ba1 - f.ba1A,
##                                     f.ea1A,   f.ea1D   =  f.ea1 - f.ea1A,
##                                     mb.a2A,   mb.a2D   =  mb.a2 - mb.a2A,
##                                     me.a2A,   me.a2D   =  me.a2 - me.a2A,
##                                     f.ba2A,   f.ba2D   =  f.ba2 - f.ba2A,
##                                     f.ea2A,   f.ea2D   =  f.ea2 - f.ea2A,
##                                     mb.A,   mb.D   =  mb.- mb.A,
##                                     me.A,   me.D   =  me. - me.A,
##                                     f.bA,   f.bD   =  f.b - f.bA,
##                                     f.eA,   f.eD   =  f.e - f.eA,
##                                     hb1b2A, hb1b2D = hb1b2 - hb1b2A, # note some of the dead cases were infected by dead partners, so can only use the index cases in any h couple
##                                     hb2b1A, hb2b1D = hb2b1 - hb2b1A,
##                                     hbeA,   hbeD   =  hbe - hbeA,
##                                     hebA,   hebD   =  heb - hebA,
##                                     hepaA,   hepaD   =  hepa - hepaA,
##                                     hpeaA,   hpeaD   =  hpea - hpeaA,
##                                     hbpaA,   hbpaD   =  hbpa - hbpaA,
##                                     hpbaA,   hpbaD   =  hpba - hpbaA,
##                                     hepA,   hepD   =  hep - hepA,
##                                     hpeA,   hpeD   =  hpe - hpeA,
##                                     hbpA,   hbpD   =  hbp - hbpA,
##                                     hpbA,   hpbD   =  hpb - hpbA,
##                                     he1e2A, he1e2D = he1e2 - he1e2A,
##                                     he2e1A, he2e1D = he2e1 - he2e1A)
##             for(ii in 1:nrow(dat)) dat$cat[ii] <- which(rmultinom(1, 1, sim.probs[ii,])==1)
##             dat$cat.nm <- names(sim.probs)[dat$cat]
##             dat$cat.nm <- factor(dat$cat.nm, levels = names(sim.probs))
##             K <- nrow(dat)
##             if(!survive) pser.a <- NA
##           }else{ ## if not simulating data calculate likelihood p(data|pars)
##           }
##       }
##     if(length(compars)==0) {
##         clprob <- NA
##         cprobs <- NA
##       }
##     if(sim) {
##         if(trace) {
##             if(give.pis) {
##                 return(list(lprob = lprob,pop.avs = pop.avs, proj12 = proj12, sim.probs, allstates = allstates,
##                             pser.a = pser.a, pser = pser, dat = dat, clprob = clprob, probs = probs, cprobs = cprobs, m.het = m.het, f.het = f.het))
##               }else{ 
##                 return(list(lprob = lprob,pop.avs = pop.avs, rrs = rrs, proj12=proj12, sim.probs, allstates = allstates,
##                             pser.a = pser.a, pser = pser, dat = dat, clprob = clprob, probs = probs, cprobs = cprobs, m.het = m.het, f.het = f.het))
##               }
##           }else{
##             return(list(lprob = lprob, pser.a = pser.a, pser = pser, dat = dat, sim.probs, allstates = allstates,
##                         clprob = clprob, probs = probs, cprobs = cprobs, m.het = m.het, f.het = f.het))
##           }
##       }else{                            # if not simulating
##         if(trace) {
##             if(give.pis) {
##                 return(list(lprob = lprob,pop.avs = pop.avs, rrs = rrs,  proj12=proj12, pis = pis, allstates = allstates,
##                             pser.a = pser.a, pser = pser, probs = probs))
##               }else{ 
##                 return(list(lprob = lprob,pop.avs = pop.avs, rrs = rrs, proj12=proj12,
##                             pser.a = pser.a, pser = pser, probs = probs))
##               }
##           }else{
##             return(list(lprob = lprob, pser.a = pser.a, pser = pser))
##           }
##       }
##   }




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
sampler <- function(sd.props = sd.props, inits,
                    acute.sc,
                    multiv = F, covar = NULL, # if multiv, sample from multivariate distribution (calculated during adaptive phase)
                    verbose = T, tell = 100, seed = 1, lrho.sd,
                    niter = 6*1000, survive,
                    nthin = 5,
                    nburn = 1000, browse=F)
  {
    if(browse)  browser()
    set.seed(seed)
    pars <- inits
    vv <- 2
    accept <- 0                  #track each parameters acceptance individually
    cur <- pcalc(pars, acute.sc = acute.sc, dat = dat, trace = T, survive = survive)     #calculate first log probability
    lprob.cur <- cur$lprob
    out <- t(data.frame(c(pars, bfp = as.numeric(pars["bmp"]*exp(pars["lrho"])),
                          cur$pop.avs, cur$rrs, cur$proj12)))
    last.it <- 0
    start <- Sys.time()
    while(vv < niter + 1)
      {
        if(verbose & vv%%tell+1==1) print(paste("on iteration",vv,"of",last.it + niter + 1))
        pars.prop <- pars              #initialize proposal parameterr vector
        ## propose new parameter vector
        if(multiv)
          {
            pars.prop <- pars.prop + rmnorm(1, mean = 0, varcov = covar)
            pars.prop <- as.vector(pars.prop) #otherwise is a matrix
            names(pars.prop) <- parnames
          }else{
            pars.prop <- pars.prop + rnorm(length(pars), mean = 0, sd = sd.props)
          }
        ## trace = T if in non-thinned iteration, or the previous one (in case of rejection)
        ## calculate proposal par log probability
        prop <- pcalc(pars.prop, acute.sc = acute.sc, dat = dat, survive = survive, trace = vv%%nthin + 1 %in% c(nthin,1))
        lprob.prop <- prop$lprob
        lmh <- lprob.prop - lprob.cur       # log Metropolis-Hastings ratio
        ## if MHR >= 1 or a uniform random # in [0,1] is <= MHR, accept otherwise reject
        if(lmh >= 0 | runif(1,0,1) <= exp(lmh))
          {
            pars <- pars.prop
            if(vv>nburn) accept <- accept + 1 #only track acceptance after burn-in
            lprob.cur <- lprob.prop
            cur <- prop
          }
        if(vv%%nthin + 1 ==1)
          {
            out <- cbind(out,t(data.frame(c(pars, bfp = as.numeric(pars["bmp"]*exp(pars["lrho"])),
                                            cur$pop.avs, cur$rrs, cur$proj12))))
          }
        vv <- vv+1
      }
    if(verbose) print(paste("took", difftime(Sys.time(),start, units = "mins"),"mins"))
    aratio <- accept/((vv-nburn))    
    return(list(out = out[,1:ncol(out)>(nburn+1)/nthin], aratio = aratio, inits = inits))
  }



## plot trajectory for selected couples
ctraj <- function(pars, dat, browse = F, width = 11, height = 8.5, last.year = 2012,
                  plot.cpls = NULL, ylab = "probability in couple serostatus / pop prevalence",
                  surv = T,             # show curves joint with survival
                  nsurv = F,            # show curves not joint with survival
                  survive = T,          # show point on surv curve
                  dead = T,          # show curve for prob of couple dying
                  cprob = T,         # show conditional probabilities
                  cols = c("orange", "darkgreen", "purple","gray"),
                  col.dead = "black",
                  lty.surv = 1, lty.dead = 1,
                  lty.nsurv = 2, blob = T, blob.cex = 3, show.susc = T, cex.leg = 1,
                  show.age = F, do.leg = T, drp = 1, x.ticks = T, t.lwd = 1,
                  do.pdf = T, pdf.name = "surv traj.pdf")
  {
    if(browse) browser()
    K <- nrow(dat)
    if(do.pdf)
      {
      pdf(pdf.name, width = width, height = height)
      par(mar=c(0,4,.5,2))
    }
    for(cpl in plot.cpls)
      {
        temp <- dat[cpl,]
        xlim <- c(80*12, (last.year-1900)*12)
        if(show.susc)
          {
            ylim <- c(-.3,1)
          }else{
            ylim <- c(-.3,.5)
          }
        plot(0,0, type = "n", xlim = xlim, axes = F,
             ylim = ylim, xlab = "",  ylab = ylab, bty = "n")
        lines(xlim[1]:xlim[2], epicm[xlim[1]:xlim[2],temp$epic.ind], col = "red", lty = 2, lwd = 3)
        lines(xlim[1]:xlim[2], epicf[xlim[1]:xlim[2],temp$epic.ind], col = "red", lty = 3, lwd = 3)
        if(surv & !nsurv & dead & do.leg)        legend("topleft", c("M+F+ & alive", "M+F- & alive","M-F+ & alive","M-F- & alive", "dead before DHS sampling"),
                                                        col = c(cols,col.dead), lwd = 3,
                                                        lty = c(rep(lty.surv,3),1,lty.dead), bty = "n", cex = cex.leg)
        if(surv & nsurv & do.leg)        legend("topleft", c("M+F+","M+F-","M-F+","M+F+ & alive",
                                                    "M+F- & alive","M-F+ & alive","M-F-"),
                                       col = c(cols[1:3],cols), lwd = 3,
                                       lty = c(rep(lty.nsurv,3),rep(lty.surv,3),1), bty = "n", cex = cex.leg)
        if(surv & !nsurv & do.leg & !dead)        legend("topleft", c("M+F+ & alive", "M+F- & alive","M-F+ & alive","M-F-"),
                                       col = cols, lwd = 3,
                                       lty = c(rep(lty.surv,3),1), bty = "n", cex = cex.leg)
        if(!surv & nsurv & do.leg)        legend("topleft", c("M+F+","M+F-","M-F+", "M-F-"),
                                       col = cols, lwd = 3,
                                       lty = c(rep(lty.nsurv,3),1), bty = "n", cex = cex.leg)
        if(do.leg) legend("left", c("HIV population prev in M", "HIV population prev in F"),
               lty = 2:3, lwd = 3, col = "red", bty = "n", cex = cex.leg)
        if(x.ticks)
          {
            axis(1, seq(xlim[1],xlim[2], by = 60), seq(xlim[1],xlim[2], by = 60)/12 + 1900, pos = 0, las = 2)
          }else{
            axis(1, seq(xlim[1],xlim[2], by = 60), labels = NA, pos = 0, las = 2)
          }
        axis(2, seq(0,1,l=5), las = 2)
        bmb <- as.numeric(pars["bmb"])
        bfb <- as.numeric(pars["bfb"])
        bme <- as.numeric(pars["bme"])
        bfe <- as.numeric(pars["bfe"])
        bmp <- as.numeric(pars["bmp"])
        rho <- exp(as.numeric(pars["lrho"]))
        bfp <- bmp * rho
        out <- data.frame(cmc = temp$tmar - temp$bd, ss = 1, mm = 0, ff = 0, hh = 0, mm.a = 0, ff.a = 0, hh.a = 0)
        ##                           pi.m.bef = 0, pi.f.bef = 0, pi.m.part = 0, pi.f.part = 0, pi.m.exc = 0, pi.f.exc = 0,
        ##                           pi.m.bef.a = 0, pi.f.bef.a = 0, pi.m.part.a = 0, pi.f.part.a = 0, pi.m.exc.a = 0, pi.f.exc.a = 0) 
        ssL <- 1
        mmL <- 0
        ffL <- 0
        hhL <- 0
        mm.aL <- 0
        ff.aL <- 0
        hh.aL <- 0
        ss <- 1
        mm <- 0
        ff <- 0
        hh <- 0
        mm.a <- 0
        ff.a <- 0
        hh.a <- 0
        for(tt in 1:max(temp$bd))
          {
            ## probabilities are non-zero only for times after started having sex and before couple formation
            m.sex <- temp$tmar-temp$bd+tt-1 >= temp$tms & temp$tmar-temp$bd+tt-1 < temp$tmar
            f.sex <- temp$tmar-temp$bd+tt-1 >= temp$tfs & temp$tmar-temp$bd+tt-1 < temp$tmar
            e.sex <- m.sex|f.sex           # either are active
            ## probability infected in month tt
            p.m.bef <- 0
            p.f.bef <- 0    
            p.m.bef.a <- 0
            p.f.bef.a <- 0    
            p.m.bef[m.sex] <- (1 - exp(-bmb * epicf[cbind(temp$tmar[m.sex]-temp$bd[m.sex]+tt-1, temp$epic.ind[m.sex])]))
            p.f.bef[f.sex] <- (1 - exp(-bfb * epicm[cbind(temp$tmar[f.sex]-temp$bd[f.sex]+tt-1, temp$epic.ind[f.sex])]))
            ## probability infected in month tt and alive at sampling
            p.m.bef.a[m.sex] <- p.m.bef[m.sex] * csurv[cbind(temp$mage[m.sex]-temp$cd[m.sex]-temp$bd[m.sex]+tt-1, temp$cd[m.sex]+temp$bd[m.sex]-tt+1)]
            p.f.bef.a[f.sex] <- p.f.bef[f.sex] * csurv[cbind(temp$fage[f.sex]-temp$cd[f.sex]-temp$bd[f.sex]+tt-1, temp$cd[f.sex]+temp$bd[f.sex]-tt+1)]
            ## iterate probabilities based on previous values for only cases where it needs uptemping
            ss[e.sex] <- ssL[e.sex]*(1-p.m.bef[e.sex])*(1-p.f.bef[e.sex])
            mm[e.sex] <- mmL[e.sex]*(1 - p.f.bef[e.sex]) + ssL[e.sex]*p.m.bef[e.sex]*(1-p.f.bef[e.sex])
            ff[e.sex] <- ffL[e.sex]*(1 - p.m.bef[e.sex]) + ssL[e.sex]*p.f.bef[e.sex]*(1-p.m.bef[e.sex])
            hh[e.sex] <- hhL[e.sex] + ssL[e.sex]*p.m.bef[e.sex]*p.f.bef[e.sex] +
              mmL[e.sex]*p.f.bef[e.sex] +
                ffL[e.sex]*p.m.bef[e.sex]
######################################################################
            ## iterate joint probabilities with survival
            mm.a[e.sex] <- mm.aL[e.sex]*(1 - p.f.bef[e.sex]) + ssL[e.sex]*p.m.bef.a[e.sex]*(1-p.f.bef[e.sex])
            ff.a[e.sex] <- ff.aL[e.sex]*(1 - p.m.bef[e.sex]) + ssL[e.sex]*p.f.bef.a[e.sex]*(1-p.m.bef[e.sex])
            hh.a[e.sex] <- hh.aL[e.sex] + ssL[e.sex]*p.m.bef.a[e.sex]*p.f.bef.a[e.sex] +
              mm.aL[e.sex]*p.f.bef.a[e.sex] +
                ff.aL[e.sex]*p.m.bef.a[e.sex]
            ## uptempe last month states
            ssL[e.sex] <- ss[e.sex]
            mmL[e.sex] <- mm[e.sex]
            ffL[e.sex] <- ff[e.sex]
            hhL[e.sex] <- hh[e.sex]
            mm.aL[e.sex] <- mm.a[e.sex]
            ff.aL[e.sex] <- ff.a[e.sex]
            hh.aL[e.sex] <- hh.a[e.sex]
            out <- rbind(out, c(cmc = temp$tmar-temp$bd+tt-1, ss, mm, ff, hh, mm.a, ff.a, hh.a))
          }
        ## probability of infection before partnership
        pi.m.bef <- mm+hh
        pi.f.bef <- ff+hh
        ## probability of infection before partnership & both individuals living to interview
        pi.m.bef.a <- mm.a+hh.a
        pi.f.bef.a <- ff.a+hh.a
        ## probability of infection by partner
        pi.m.part <- 0
        pi.f.part <- 0
        ## probability of infection by partner & both individuals living to interview
        pi.m.part.a <- 0
        pi.f.part.a <- 0
        ## probability of infection extracouply
        pi.m.exc <- 0
        pi.f.exc <- 0
        ## probability of infection extracouply & both individuals living to interview
        pi.m.exc.a <- 0
        pi.f.exc.a <- 0
        ## Track serodiscordant couples by route of infection for bookkeeping later
        mm.bef.a <- mm.a
        ff.bef.a <- ff.a    
        mm.exc.a <- 0
        ff.exc.a <- 0
        mm.bef.aL <- mm.a
        ff.bef.aL <- ff.a    
        mm.exc.aL <- 0
        ff.exc.aL <- 0
        ## probability of being infected by partner (constant, used inside loop)
        p.m.part <- 1 - exp(-bmp)
        p.f.part <- 1 - exp(-bfp)
        ## Now loop through marriage
        for(tt in 1:max(temp$cd-1))
          {
            ## are partners formed in a couple?
            fmd <- temp$cd >= tt              
######################################################################
            ## everything below is automatically sum(fmd) length except p.m/f.part which are length 1
            ## Survival probabilities
            s.p.m <- csurv[cbind(temp$mage[fmd]-temp$cd[fmd]+tt-1, temp$tint[fmd] - (temp$tmar[fmd] + tt - 1))]
            s.p.f <- csurv[cbind(temp$fage[fmd]-temp$cd[fmd]+tt-1, temp$tint[fmd] - (temp$tmar[fmd] + tt - 1))]
            ## Transmission probabilities from partner (jointly with survival)
            p.m.part.a <- p.m.part * s.p.m
            p.f.part.a <- p.f.part * s.p.f
            ## probability infected extracouply in the ttc-th month of couple
            p.m.exc <- (1 - exp(-bme*epicf[cbind(temp$tmar[fmd]+(tt-1), temp$epic.ind[fmd])]))
            p.f.exc <- (1 - exp(-bfe*epicm[cbind(temp$tmar[fmd]+(tt-1), temp$epic.ind[fmd])]))
            p.m.exc.a <- p.m.exc * s.p.m
            p.f.exc.a <- p.f.exc * s.p.f
######################################################################
            ## iterate probabilities
            ss[fmd] <- ssL[fmd]*(1-p.m.exc)*(1-p.f.exc)
            mm[fmd] <- mmL[fmd]*(1-p.f.exc)*(1-p.f.part) + ssL[fmd]*p.m.exc*(1-p.f.exc)
            ff[fmd] <- ffL[fmd]*(1-p.m.exc)*(1-p.m.part) + ssL[fmd]*p.f.exc*(1-p.m.exc)
            hh[fmd] <- hhL[fmd] + ssL[fmd]* p.m.exc*p.f.exc +
              mmL[fmd]*(p.f.part + (1-p.f.part)*p.f.exc) +
                ffL[fmd]*(p.m.part + (1-p.m.part)*p.m.exc)
######################################################################
            ## iterate probabilities jointly with survival until survey
            ## Note for probabilities of not being infected, we don't use the joint probability with being alive at sampling
            mm.a[fmd] <- mm.aL[fmd]*(1-p.f.exc)*(1-p.f.part) + ssL[fmd]*p.m.exc.a*(1-p.f.exc)
            ff.a[fmd] <- ff.aL[fmd]*(1-p.m.exc)*(1-p.m.part) + ssL[fmd]*p.f.exc.a*(1-p.m.exc)
            hh.a[fmd] <- hh.aL[fmd] + ssL[fmd]  *p.m.exc.a*p.f.exc.a +
              mm.aL[fmd]*(p.f.part.a + (1-p.f.part)*p.f.exc.a) +
                ff.aL[fmd]*(p.m.part.a + (1-p.m.part)*p.m.exc.a)
            ## Track cumulative incidence probabilities for partner infections
            pi.m.part[fmd] <- pi.m.part[fmd] + ffL[fmd]*p.m.part
            pi.f.part[fmd] <- pi.f.part[fmd] + mmL[fmd]*p.f.part        
            pi.m.part.a[fmd] <- pi.m.part.a[fmd] + ff.aL[fmd]*p.m.part.a
            pi.f.part.a[fmd] <- pi.f.part.a[fmd] + mm.aL[fmd]*p.f.part.a
            ## Track how many serodiscordant couples there are with
            ## infections from before partnership 
            mm.bef.a[fmd] <- mm.bef.aL[fmd]*(1-p.f.exc)*(1-p.f.part)
            ff.bef.a[fmd] <- ff.bef.aL[fmd]*(1-p.m.exc)*(1-p.m.part)
            ## Track cumulative incidence probabilities for infections
            ## prior to partnership just subtracting off partners with
            ## incident infections that die and remove couple from
            ## sample-able couples
            pi.m.bef.a[fmd] <- pi.m.bef.a[fmd] - mm.bef.aL[fmd]*(p.f.part + (1-p.f.part)*p.f.exc)*(1-s.p.f)
            pi.f.bef.a[fmd] <- pi.f.bef.a[fmd] - ff.bef.aL[fmd]*(p.m.part + (1-p.m.part)*p.m.exc)*(1-s.p.m)
            ## Track cumulative incidence probabilties for extracouple infections
            mm.exc.a[fmd] <- mm.aL[fmd]-mm.bef.aL[fmd]
            ff.exc.a[fmd] <- ff.aL[fmd]-ff.bef.aL[fmd]        
            pi.m.exc[fmd] <- pi.m.exc[fmd] + (ssL[fmd] + ffL[fmd]*(1-p.m.part))*p.m.exc
            pi.f.exc[fmd] <- pi.f.exc[fmd] + (ssL[fmd] + mmL[fmd]*(1-p.f.part))*p.f.exc
            pi.m.exc.a[fmd] <- pi.m.exc.a[fmd] + (ssL[fmd] + ff.aL[fmd]*(1-p.m.part))*p.m.exc.a -
              mm.exc.aL[fmd]*(1-p.f.part)*p.f.exc*(1-s.p.f)
            pi.f.exc.a[fmd] <- pi.f.exc.a[fmd] + (ssL[fmd] + mm.aL[fmd]*(1-p.f.part))*p.f.exc.a -
              ff.exc.aL[fmd]*(1-p.m.part)*p.m.exc*(1-s.p.m)       
            ## update last month states
            ssL[fmd] <- ss[fmd]
            mmL[fmd] <- mm[fmd]
            ffL[fmd] <- ff[fmd]
            hhL[fmd] <- hh[fmd]       
            mm.aL[fmd] <- mm.a[fmd]
            ff.aL[fmd] <- ff.a[fmd]
            hh.aL[fmd] <- hh.a[fmd]
            mm.exc.aL[fmd] <- mm.exc.a[fmd]
            ff.exc.aL[fmd] <- ff.exc.a[fmd]            
            mm.bef.aL[fmd] <- mm.bef.a[fmd]
            ff.bef.aL[fmd] <- ff.bef.a[fmd]
            out <- rbind(out, c(cmc = temp$tmar+tt-1, ss, mm, ff, hh, mm.a, ff.a, hh.a))
          }
        pser.a <- cbind(hh.a, mm.a, ff.a, ss)
        pser <- cbind(hh, mm, ff, ss)
        if(cprob)
          {
            scalar <- out$mm.a + out$ff.a + out$hh.a + out$ss
          }else{
            scalar <- 1
          }
        if(show.susc)   lines(out$cmc, out$ss/scalar, col = cols[4], lty = 1, lwd = 2)
        if(nsurv)
          {
            lines(out$cmc, out$mm, col = cols[2], lty = lty.nsurv, lwd = 2)
            lines(out$cmc, out$ff, col = cols[3], lty = lty.nsurv, lwd = 2)
            lines(out$cmc, out$hh, col = cols[1], lty = lty.nsurv, lwd = 2)
            if(!survive & blob) points(out$cmc[nrow(out)], pser[temp$ser], pch = 19,
                               col = cols[temp$ser], cex = 3)
          }
        if(surv)
          {
            lines(out$cmc, out$mm.a/scalar, col = cols[2], lty = lty.surv, lwd = 2)
            lines(out$cmc, out$ff.a/scalar, col = cols[3], lty = lty.surv, lwd = 2)
            lines(out$cmc, out$hh.a/scalar, col = cols[1], lty = lty.surv, lwd = 2)
            if(survive & blob) points(out$cmc[nrow(out)], pser.a[temp$ser]/scalar[nrow(out)], pch = 19,
                   col = cols[temp$ser], cex = blob.cex)
          }
        if(dead)
          {
             lines(out$cmc, out$mm + out$ff + out$hh - out$mm.a - out$ff.a - out$hh.a, col = col.dead, lty = lty.dead, lwd = 2)
           }
        segments(temp$tmar,-.01,temp$tmar,-drp*.1, lwd = t.lwd)
        text(temp$tmar,-drp*.1,"couple formation",pos=1)
        segments(temp$tms,-.01,temp$tms,-drp*.2, col = cols[2], lwd = t.lwd)
        text(temp$tms,-drp*.2,"male first sex",pos=2, col = cols[2])
        segments(temp$tfs,-.01,temp$tfs,-drp*.3, col = cols[3], lwd = t.lwd)
        text(temp$tfs,-drp*.3,"female first sex",pos=2, col = cols[3])
        segments(temp$tint,-.01,temp$tint,-drp*.2, lwd = t.lwd)
        text(temp$tint,-drp*.2,"DHS survey",pos=4, col = "black")
        if(show.age)
          {
            mtext(paste("male is", round(temp$mage/12,0),"yrs"), adj = .2, side = 3, line = -3)
            mtext(paste("female is", round(temp$fage/12,0),"yrs"), adj = .2, side = 3, line = -5)
          }
      }
    if(do.pdf) dev.off()
  }


## plot trajectory for selected couples
ctraj.area <- function(pars, dat, browse = F, width = 11, height = 8.5, last.year = 2012,
                       plot.cpls = NULL, ylab = "probability in couple serostatus / pop prevalence",
                       surv = T,             # show curves joint with survival
                       nsurv = F,            # show curves not joint with survival
                       survive = T,          # show point on surv curve
                       dead = T,          # show curve for prob of couple dying
                       cprob = T,         # show conditional probabilities
                       cols = c("orange", "darkgreen", "purple","gray"), cex.t = .8,
                       col.dead = "black",
                       lty.surv = 1, lty.dead = 1,
                       lty.nsurv = 2, blob = T, blob.cex = 3, show.susc = T, cex.leg = 1,
                       show.age = F, do.leg = T, drp = 1, x.ticks = T, t.lwd = 1,
                       do.pdf = T, pdf.name = "surv traj.pdf")
  {
    if(browse) browser()
    K <- nrow(dat)
    if(do.pdf)
      {
        pdf(pdf.name, width = width, height = height)
        par(mar=c(0,4,.5,2))
      }
    for(cpl in plot.cpls)
      {
        temp <- dat[cpl,]
        xlim <- c(80*12, (last.year-1900)*12)
        if(show.susc)
          {
            ylim <- c(-.3,1)
          }else{
            ylim <- c(-.3,.5)
          }
        plot(0,0, type = "n", xlim = xlim, axes = F,
             ylim = ylim, xlab = "",  ylab = ylab, bty = "n")
        if(surv & !nsurv & dead & do.leg)        legend("topleft", c("M+F+ & alive", "M+F- & alive","M-F+ & alive","M-F- & alive", "dead before DHS sampling"),
                                                        col = c(cols,col.dead), lwd = 3,
                                                        lty = c(rep(lty.surv,3),1,lty.dead), bty = "n", cex = cex.leg)
        if(surv & nsurv & do.leg)        legend("topleft", c("M+F+","M+F-","M-F+","M+F+ & alive",
                                                             "M+F- & alive","M-F+ & alive","M-F-"),
                                                col = c(cols[1:3],cols), lwd = 3,
                                                lty = c(rep(lty.nsurv,3),rep(lty.surv,3),1), bty = "n", cex = cex.leg)
        if(surv & !nsurv & do.leg & !dead)        legend("topleft", c("M+F+ & alive", "M+F- & alive","M-F+ & alive","M-F-"),
                                                         col = cols, lwd = 3,
                                                         lty = c(rep(lty.surv,3),1), bty = "n", cex = cex.leg)
        if(!surv & nsurv & do.leg)        legend("topleft", c("M+F+","M+F-","M-F+", "M-F-"),
                                                 col = cols, lwd = 3,
                                                 lty = c(rep(lty.nsurv,3),1), bty = "n", cex = cex.leg)
        if(do.leg) legend("left", c("HIV population prev in M", "HIV population prev in F"),
                          lty = 2:3, lwd = 3, col = "red", bty = "n", cex = cex.leg)
        if(x.ticks)
          {
            axis(1, seq(xlim[1],xlim[2], by = 12*10), seq(xlim[1],xlim[2], by = 12*10)/12 + 1900, pos = 0, las = 2, cex.axis = cex.leg)
            axis(1, seq(xlim[1],xlim[2], by = 60), labels = NA, pos = 0, las = 2)            
          }else{
            axis(1, seq(xlim[1],xlim[2], by = 60), labels = NA, pos = 0, las = 2)
          }
        axis(2, seq(0,1,l=3), las = 2, cex.axis = cex.leg)
        axis(2, seq(0,1,l=5), labels = NA)
        bmb <- as.numeric(pars["bmb"])
        bfb <- as.numeric(pars["bfb"])
        bme <- as.numeric(pars["bme"])
        bfe <- as.numeric(pars["bfe"])
        bmp <- as.numeric(pars["bmp"])
        rho <- exp(as.numeric(pars["lrho"]))
        bfp <- bmp * rho
        out <- data.frame(cmc = temp$tmar - temp$bd, ss = 1, mm = 0, ff = 0, hh = 0, mm.a = 0, ff.a = 0, hh.a = 0)
        ##                           pi.m.bef = 0, pi.f.bef = 0, pi.m.part = 0, pi.f.part = 0, pi.m.exc = 0, pi.f.exc = 0,
        ##                           pi.m.bef.a = 0, pi.f.bef.a = 0, pi.m.part.a = 0, pi.f.part.a = 0, pi.m.exc.a = 0, pi.f.exc.a = 0) 
        ssL <- 1
        mmL <- 0
        ffL <- 0
        hhL <- 0
        mm.aL <- 0
        ff.aL <- 0
        hh.aL <- 0
        ss <- 1
        mm <- 0
        ff <- 0
        hh <- 0
        mm.a <- 0
        ff.a <- 0
        hh.a <- 0
        for(tt in 1:max(temp$bd))
          {
            ## probabilities are non-zero only for times after started having sex and before couple formation
            m.sex <- temp$tmar-temp$bd+tt-1 >= temp$tms & temp$tmar-temp$bd+tt-1 < temp$tmar
            f.sex <- temp$tmar-temp$bd+tt-1 >= temp$tfs & temp$tmar-temp$bd+tt-1 < temp$tmar
            e.sex <- m.sex|f.sex           # either are active
            ## probability infected in month tt
            p.m.bef <- 0
            p.f.bef <- 0    
            p.m.bef.a <- 0
            p.f.bef.a <- 0    
            p.m.bef[m.sex] <- (1 - exp(-bmb * epicf[cbind(temp$tmar[m.sex]-temp$bd[m.sex]+tt-1, temp$epic.ind[m.sex])]))
            p.f.bef[f.sex] <- (1 - exp(-bfb * epicm[cbind(temp$tmar[f.sex]-temp$bd[f.sex]+tt-1, temp$epic.ind[f.sex])]))
            ## probability infected in month tt and alive at sampling
            p.m.bef.a[m.sex] <- p.m.bef[m.sex] * csurv[cbind(temp$mage[m.sex]-temp$cd[m.sex]-temp$bd[m.sex]+tt-1, temp$cd[m.sex]+temp$bd[m.sex]-tt+1)]
            p.f.bef.a[f.sex] <- p.f.bef[f.sex] * csurv[cbind(temp$fage[f.sex]-temp$cd[f.sex]-temp$bd[f.sex]+tt-1, temp$cd[f.sex]+temp$bd[f.sex]-tt+1)]
            ## iterate probabilities based on previous values for only cases where it needs uptemping
            ss[e.sex] <- ssL[e.sex]*(1-p.m.bef[e.sex])*(1-p.f.bef[e.sex])
            mm[e.sex] <- mmL[e.sex]*(1 - p.f.bef[e.sex]) + ssL[e.sex]*p.m.bef[e.sex]*(1-p.f.bef[e.sex])
            ff[e.sex] <- ffL[e.sex]*(1 - p.m.bef[e.sex]) + ssL[e.sex]*p.f.bef[e.sex]*(1-p.m.bef[e.sex])
            hh[e.sex] <- hhL[e.sex] + ssL[e.sex]*p.m.bef[e.sex]*p.f.bef[e.sex] +
              mmL[e.sex]*p.f.bef[e.sex] +
                ffL[e.sex]*p.m.bef[e.sex]
######################################################################
            ## iterate joint probabilities with survival
            mm.a[e.sex] <- mm.aL[e.sex]*(1 - p.f.bef[e.sex]) + ssL[e.sex]*p.m.bef.a[e.sex]*(1-p.f.bef[e.sex])
            ff.a[e.sex] <- ff.aL[e.sex]*(1 - p.m.bef[e.sex]) + ssL[e.sex]*p.f.bef.a[e.sex]*(1-p.m.bef[e.sex])
            hh.a[e.sex] <- hh.aL[e.sex] + ssL[e.sex]*p.m.bef.a[e.sex]*p.f.bef.a[e.sex] +
              mm.aL[e.sex]*p.f.bef.a[e.sex] +
                ff.aL[e.sex]*p.m.bef.a[e.sex]
            ## uptempe last month states
            ssL[e.sex] <- ss[e.sex]
            mmL[e.sex] <- mm[e.sex]
            ffL[e.sex] <- ff[e.sex]
            hhL[e.sex] <- hh[e.sex]
            mm.aL[e.sex] <- mm.a[e.sex]
            ff.aL[e.sex] <- ff.a[e.sex]
            hh.aL[e.sex] <- hh.a[e.sex]
            out <- rbind(out, c(cmc = temp$tmar-temp$bd+tt-1, ss, mm, ff, hh, mm.a, ff.a, hh.a))
          }
        ## probability of infection before partnership
        pi.m.bef <- mm+hh
        pi.f.bef <- ff+hh
        ## probability of infection before partnership & both individuals living to interview
        pi.m.bef.a <- mm.a+hh.a
        pi.f.bef.a <- ff.a+hh.a
        ## probability of infection by partner
        pi.m.part <- 0
        pi.f.part <- 0
        ## probability of infection by partner & both individuals living to interview
        pi.m.part.a <- 0
        pi.f.part.a <- 0
        ## probability of infection extracouply
        pi.m.exc <- 0
        pi.f.exc <- 0
        ## probability of infection extracouply & both individuals living to interview
        pi.m.exc.a <- 0
        pi.f.exc.a <- 0
        ## Track serodiscordant couples by route of infection for bookkeeping later
        mm.bef.a <- mm.a
        ff.bef.a <- ff.a    
        mm.exc.a <- 0
        ff.exc.a <- 0
        mm.bef.aL <- mm.a
        ff.bef.aL <- ff.a    
        mm.exc.aL <- 0
        ff.exc.aL <- 0
        ## probability of being infected by partner (constant, used inside loop)
        p.m.part <- 1 - exp(-bmp)
        p.f.part <- 1 - exp(-bfp)
        ## Now loop through marriage
        for(tt in 1:max(temp$cd-1))
          {
            ## are partners formed in a couple?
            fmd <- temp$cd >= tt              
######################################################################
            ## everything below is automatically sum(fmd) length except p.m/f.part which are length 1
            ## Survival probabilities
            s.p.m <- csurv[cbind(temp$mage[fmd]-temp$cd[fmd]+tt-1, temp$tint[fmd] - (temp$tmar[fmd] + tt - 1))]
            s.p.f <- csurv[cbind(temp$fage[fmd]-temp$cd[fmd]+tt-1, temp$tint[fmd] - (temp$tmar[fmd] + tt - 1))]
            ## Transmission probabilities from partner (jointly with survival)
            p.m.part.a <- p.m.part * s.p.m
            p.f.part.a <- p.f.part * s.p.f
            ## probability infected extracouply in the ttc-th month of couple
            p.m.exc <- (1 - exp(-bme*epicf[cbind(temp$tmar[fmd]+(tt-1), temp$epic.ind[fmd])]))
            p.f.exc <- (1 - exp(-bfe*epicm[cbind(temp$tmar[fmd]+(tt-1), temp$epic.ind[fmd])]))
            p.m.exc.a <- p.m.exc * s.p.m
            p.f.exc.a <- p.f.exc * s.p.f
######################################################################
            ## iterate probabilities
            ss[fmd] <- ssL[fmd]*(1-p.m.exc)*(1-p.f.exc)
            mm[fmd] <- mmL[fmd]*(1-p.f.exc)*(1-p.f.part) + ssL[fmd]*p.m.exc*(1-p.f.exc)
            ff[fmd] <- ffL[fmd]*(1-p.m.exc)*(1-p.m.part) + ssL[fmd]*p.f.exc*(1-p.m.exc)
            hh[fmd] <- hhL[fmd] + ssL[fmd]* p.m.exc*p.f.exc +
              mmL[fmd]*(p.f.part + (1-p.f.part)*p.f.exc) +
                ffL[fmd]*(p.m.part + (1-p.m.part)*p.m.exc)
######################################################################
            ## iterate probabilities jointly with survival until survey
            ## Note for probabilities of not being infected, we don't use the joint probability with being alive at sampling
            mm.a[fmd] <- mm.aL[fmd]*(1-p.f.exc)*(1-p.f.part) + ssL[fmd]*p.m.exc.a*(1-p.f.exc)
            ff.a[fmd] <- ff.aL[fmd]*(1-p.m.exc)*(1-p.m.part) + ssL[fmd]*p.f.exc.a*(1-p.m.exc)
            hh.a[fmd] <- hh.aL[fmd] + ssL[fmd]  *p.m.exc.a*p.f.exc.a +
              mm.aL[fmd]*(p.f.part.a + (1-p.f.part)*p.f.exc.a) +
                ff.aL[fmd]*(p.m.part.a + (1-p.m.part)*p.m.exc.a)
            ## Track cumulative incidence probabilities for partner infections
            pi.m.part[fmd] <- pi.m.part[fmd] + ffL[fmd]*p.m.part
            pi.f.part[fmd] <- pi.f.part[fmd] + mmL[fmd]*p.f.part        
            pi.m.part.a[fmd] <- pi.m.part.a[fmd] + ff.aL[fmd]*p.m.part.a
            pi.f.part.a[fmd] <- pi.f.part.a[fmd] + mm.aL[fmd]*p.f.part.a
            ## Track how many serodiscordant couples there are with
            ## infections from before partnership 
            mm.bef.a[fmd] <- mm.bef.aL[fmd]*(1-p.f.exc)*(1-p.f.part)
            ff.bef.a[fmd] <- ff.bef.aL[fmd]*(1-p.m.exc)*(1-p.m.part)
            ## Track cumulative incidence probabilities for infections
            ## prior to partnership just subtracting off partners with
            ## incident infections that die and remove couple from
            ## sample-able couples
            pi.m.bef.a[fmd] <- pi.m.bef.a[fmd] - mm.bef.aL[fmd]*(p.f.part + (1-p.f.part)*p.f.exc)*(1-s.p.f)
            pi.f.bef.a[fmd] <- pi.f.bef.a[fmd] - ff.bef.aL[fmd]*(p.m.part + (1-p.m.part)*p.m.exc)*(1-s.p.m)
            ## Track cumulative incidence probabilties for extracouple infections
            mm.exc.a[fmd] <- mm.aL[fmd]-mm.bef.aL[fmd]
            ff.exc.a[fmd] <- ff.aL[fmd]-ff.bef.aL[fmd]        
            pi.m.exc[fmd] <- pi.m.exc[fmd] + (ssL[fmd] + ffL[fmd]*(1-p.m.part))*p.m.exc
            pi.f.exc[fmd] <- pi.f.exc[fmd] + (ssL[fmd] + mmL[fmd]*(1-p.f.part))*p.f.exc
            pi.m.exc.a[fmd] <- pi.m.exc.a[fmd] + (ssL[fmd] + ff.aL[fmd]*(1-p.m.part))*p.m.exc.a -
              mm.exc.aL[fmd]*(1-p.f.part)*p.f.exc*(1-s.p.f)
            pi.f.exc.a[fmd] <- pi.f.exc.a[fmd] + (ssL[fmd] + mm.aL[fmd]*(1-p.f.part))*p.f.exc.a -
              ff.exc.aL[fmd]*(1-p.m.part)*p.m.exc*(1-s.p.m)       
            ## update last month states
            ssL[fmd] <- ss[fmd]
            mmL[fmd] <- mm[fmd]
            ffL[fmd] <- ff[fmd]
            hhL[fmd] <- hh[fmd]       
            mm.aL[fmd] <- mm.a[fmd]
            ff.aL[fmd] <- ff.a[fmd]
            hh.aL[fmd] <- hh.a[fmd]
            mm.exc.aL[fmd] <- mm.exc.a[fmd]
            ff.exc.aL[fmd] <- ff.exc.a[fmd]            
            mm.bef.aL[fmd] <- mm.bef.a[fmd]
            ff.bef.aL[fmd] <- ff.bef.a[fmd]
            out <- rbind(out, c(cmc = temp$tmar+tt-1, ss, mm, ff, hh, mm.a, ff.a, hh.a))
          }
        pser.a <- cbind(hh.a, mm.a, ff.a, ss)
        pser <- cbind(hh, mm, ff, ss)
        if(cprob)
          {
            scalar <- out$mm.a + out$ff.a + out$hh.a + out$ss
          }else{
            scalar <- 1
          }
                                        #        if(show.susc)   lines(out$cmc, out$ss/scalar, col = cols[4], lty = 1, lwd = 2)
        if(nsurv)
          {
            lines(out$cmc, out$mm, col = cols[2], lty = lty.nsurv, lwd = 2)
            lines(out$cmc, out$ff, col = cols[3], lty = lty.nsurv, lwd = 2)
            lines(out$cmc, out$hh, col = cols[1], lty = lty.nsurv, lwd = 2)
            if(!survive & blob) points(out$cmc[nrow(out)], pser[temp$ser], pch = 19,
                                       col = cols[temp$ser], cex = 3)
          }
        if(surv)
          {
            if(out$cmc[1] > xlim[1])
              {
                out.begin <- out[rep(1, out$cmc[1]-xlim[1]),]
                out.begin$cmc <- xlim[1]:(out$cmc[1]-1)
                outpl <- rbind(out.begin,out)
              }else{
                outpl <- out[out$cmc >= xlim[1],]
              }
            scalar <- outpl$mm.a + outpl$ff.a + outpl$hh.a + outpl$ss
            polygon(c(outpl$cmc, rev(outpl$cmc)), c(outpl$hh.a/scalar, rep(0, nrow(outpl))), col = cols[1])
            polygon(c(outpl$cmc, rev(outpl$cmc)), c(outpl$hh.a/scalar, rev(outpl$hh.a/scalar + outpl$ff.a/scalar)), col = cols[3])
            polygon(c(outpl$cmc, rev(outpl$cmc)), c(outpl$hh.a/scalar + outpl$ff.a/scalar, rev(outpl$hh.a/scalar + outpl$ff.a/scalar + outpl$mm.a/scalar)),  col = cols[2])
            if(show.susc)   polygon(c(outpl$cmc, rev(outpl$cmc)), c(outpl$hh.a/scalar + outpl$ff.a/scalar + outpl$mm.a/scalar,
                                                                rev(outpl$hh.a/scalar + outpl$ff.a/scalar + outpl$mm.a/scalar + outpl$ss/scalar)), col = cols[4])
            if(survive & blob) points(outpl$cmc[nrow(outpl)]+12*2, pser.a[temp$ser]/scalar[nrow(outpl)], pch = 19, lwd = 1,
                   col = cols[temp$ser], cex = blob.cex)
            if(survive & blob) points(outpl$cmc[nrow(outpl)]+12*2, pser.a[temp$ser]/scalar[nrow(outpl)], pch = 21, lwd = 1,
                   col = "black", cex = blob.cex)
          }
        if(dead)
          {
             lines(outpl$cmc, outpl$mm + outpl$ff + outpl$hh - outpl$mm.a - outpl$ff.a - outpl$hh.a, col = col.dead, lty = lty.dead, lwd = 2)
           }
        segments(temp$tmar,-.01,temp$tmar,-drp, lwd = t.lwd)
        text(temp$tmar,-drp,expression(t),pos=1, cex = cex.t)
        segments(temp$tms,-.01,temp$tms,-drp, col = cols[2], lwd = t.lwd)
        text(temp$tms,-drp,expression(t),pos=1, col = cols[2], cex = cex.t)
        segments(temp$tfs,-.01,temp$tfs,-drp, col = cols[3], lwd = t.lwd)
        text(temp$tfs,-drp,expression(t),pos=1, col = cols[3], cex = cex.t)
        segments(temp$tint,-.01,temp$tint,-drp, lwd = t.lwd)
        text(temp$tint,-drp,expression(t),pos=1, col = "black", cex = cex.t)
#        lines(xlim[1]:xlim[2], epicm[xlim[1]:xlim[2],temp$epic.ind], col = "red", lty = 2, lwd = 3)
 #       lines(xlim[1]:xlim[2], epicf[xlim[1]:xlim[2],temp$epic.ind], col = "red", lty = 3, lwd = 3)
        if(show.age)
          {
            mtext(paste("male is", round(temp$mage/12,0),"yrs"), adj = .2, side = 3, line = -3)
            mtext(paste("female is", round(temp$fage/12,0),"yrs"), adj = .2, side = 3, line = -5)
          }
      }
    if(do.pdf) dev.off()
  }

## Plot posterior pairwise-correlations & histograms
sbpairs <- function(posts, truepars = NULL, file.nm, width = 10, height = 10, show.lines = T,
                    cex = 1, col = "black", nrpoints = 200, do.pdf = F, do.jpeg = T, big.ranges = T,
                    cex.axis = 1.5, cex.nm = 2.5, greek = T, show.cis = T, browse = F)
  {
    if(browse) browser()
    if(do.pdf) pdf(paste(file.nm, ".pdf", sep=""), width = width, height = height)
    if(do.jpeg) jpeg(paste(file.nm, ".jpeg", sep=""), width = width*100, height = height*100)
    ranges <- apply(rbind(posts,truepars), 2, range)
    if(big.ranges)
      {
        ranges[1,1:5] <- 0
        ranges[2,1:5] <- .06 #1.5* ranges[2,1:5]
        ranges[1,6] <- ranges[1,6] - .5
        ranges[2,6] <- ranges[2,6] + .5
      }
    cis <- apply(rbind(posts,truepars), 2, function(x) quantile(x, c(.025, .975)))
    par(mar=rep(3,4),oma=rep(2,4),mfrow=rep(ncol(posts),2))
    parnames <- c(expression(beta[Mb]), expression(beta[Fb]),
                  expression(beta[Me]), expression(beta[Fe]),
                  expression(beta[Mp]), expression(log(rho)))
    for(ii in 1:ncol(posts))
      {
        for(jj in 1:ncol(posts))
          {
            if(ii==jj)
              {
                hist(posts[,ii], breaks = 40, col = "black", main = "", xlab = "", ylab = "",
                     xlim = ranges[,ii], las = 2, cex.axis = cex.axis)
                if(show.lines) abline(v=truepars[ii], col = "red", lwd = 3)
                if(show.cis) abline(v=cis[,ii], col = "yellow", lwd = 3)
                mtext(parnames[ii], side = 3, line = -2, adj = .98, col = "red", cex = cex.nm)
              }else
            {
              smoothScatter(posts[,jj],posts[,ii], cex = cex, col = col, las = 2, cex.axis = cex.axis,
                            main = "", xlab = "", ylab = "", nrpoints = nrpoints,
                            xlim = ranges[,jj], ylim = ranges[,ii])
              if(show.lines) abline(v=truepars[jj], col = "red", lwd = 2)
              if(show.lines) abline(h=truepars[ii], col = "red", lwd = 2)              
            }
          }
      }
    if(do.pdf | do.jpeg) dev.off()
  }


## Make scatterplot with marginal histograms
    scatterhist = function(x, y, xlab="", ylab="", xb = 5, yb = .1, border = NA, cex = .4, xlim = c(0,60), ylim = c(0,1), ..., nn=F,
      col = "black", browse = F)
      {
        if(browse) browser()
        par(oma=c(1,1,0,0))        
        zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
        layout(zones, widths=c(3.5,2), heights=c(2,3.5))
        xbrks <- seq(min(xlim),max(xlim), by = xb)
        ybrks <- seq(min(ylim),max(ylim), by = yb)        
        xhist = hist(x, plot=FALSE, breaks = xbrks)
        yhist = hist(y, plot=FALSE, breaks = ybrks)
        top = max(c(xhist$counts, yhist$counts))
        if(nn) top <- max(xhist$counts)
        par(mar=c(2,2,1,1))
        if(!nn)
          {            
            plot(x,y, bty = "n", ylim = ylim, xlim = xlim, pch = 19, cex = cex, las = 2, axes=T, col = col)
          }else{
            plot(x,y, bty = "n", axes = F, type = "n")
          }
        par(mar=c(0,2,1,1))
        bp <- barplot(xhist$counts, axes=F, ylim=c(0, top), space=0, border = border, col = "black")
        show <- 1:length(xbrks) %% 2 == 1
        axis(1, (c(bp,max(bp)+1)-.5)[show], NA, las = 2)
        if(!nn)
          {
            par(mar=c(2,0,1,1))
            bp <- barplot(yhist$counts, axes=F, xlim=c(0, top), space=0, horiz=TRUE, border = border, col = "black")
            show <- 1:length(ybrks) %% 2 == 1
            axis(2, (c(bp,max(bp)+1)-.5)[show], NA, las = 2)
          }
      }


    add.cat <- function(dat)
      { # add categories corresponding to routes of infection for both dead & alive couples
        cat <- rep(NA, nrow(dat))
        cat[dat$ser==4] <- 's..'
        cat[dat$ser==2 & dat$mdoi < dat$tmar] <- 'mb.'
        cat[dat$ser==2 & dat$mdoi >= dat$tmar] <- 'me.'
        cat[dat$ser==3 & dat$fdoi < dat$tmar] <- 'f.b'
        cat[dat$ser==3 & dat$fdoi >= dat$tmar] <- 'f.e'
        cat[dat$ser==1 & dat$mdoi <= dat$fdoi & dat$mcoi=='b' & dat$fcoi=='b'] <- 'hb1b2'
        cat[dat$ser==1 & dat$mdoi > dat$fdoi & dat$mcoi=='b' & dat$fcoi=='b'] <- 'hb2b1'
        cat[dat$ser==1 & dat$mdoi <= dat$fdoi & dat$mcoi=='e' & dat$fcoi=='e'] <- 'he1e2'
        cat[dat$ser==1 & dat$mdoi > dat$fdoi & dat$mcoi=='e' & dat$fcoi=='e'] <- 'he2e1'
        cat[dat$ser==1 & dat$mcoi=='b' & dat$fcoi=='e'] <- 'hbe'
        cat[dat$ser==1 & dat$mcoi=='e' & dat$fcoi=='b'] <- 'heb'
        cat[dat$ser==1 & dat$mcoi=='b' & dat$fcoi=='p' & dat$fcoi.phase=='a'] <- 'hbpa'
        cat[dat$ser==1 & dat$mcoi=='b' & dat$fcoi=='p' & !dat$fcoi.phase=='a'] <- 'hbp'
        cat[dat$ser==1 & dat$mcoi=='p' & dat$fcoi=='b' & dat$mcoi.phase=='a'] <- 'hpba'
        cat[dat$ser==1 & dat$mcoi=='p' & dat$fcoi=='b' & !dat$mcoi.phase=='a'] <- 'hpb'
        cat[dat$ser==1 & dat$mcoi=='e' & dat$fcoi=='p' & dat$fcoi.phase=='a'] <- 'hepa'
        cat[dat$ser==1 & dat$mcoi=='e' & dat$fcoi=='p' & !dat$fcoi.phase=='a'] <- 'hep'
        cat[dat$ser==1 & dat$mcoi=='p' & dat$fcoi=='e' & dat$mcoi.phase=='a'] <- 'hpea'
        cat[dat$ser==1 & dat$mcoi=='p' & dat$fcoi=='e' & !dat$mcoi.phase=='a'] <- 'hpe'
        ## and survived
        catA <- rep(NA, nrow(dat))
        catA[dat$alive & dat$ser==4] <- 's..'
        catA[dat$alive & dat$ser==2 & dat$mdoi < dat$tmar] <- 'mb.A'
        catA[dat$alive & dat$ser==2 & dat$mdoi >= dat$tmar] <- 'me.A'
        catA[dat$alive & dat$ser==3 & dat$fdoi < dat$tmar] <- 'f.bA'
        catA[dat$alive & dat$ser==3 & dat$fdoi >= dat$tmar] <- 'f.eA'
        catA[dat$alive & dat$ser==1 & dat$mdoi <= dat$fdoi & dat$mcoi=='b' & dat$fcoi=='b'] <- 'hb1b2A'
        catA[dat$alive & dat$ser==1 & dat$mdoi > dat$fdoi & dat$mcoi=='b' & dat$fcoi=='b'] <- 'hb2b1A'
        catA[dat$alive & dat$ser==1 & dat$mdoi <= dat$fdoi & dat$mcoi=='e' & dat$fcoi=='e'] <- 'he1e2A'
        catA[dat$alive & dat$ser==1 & dat$mdoi > dat$fdoi & dat$mcoi=='e' & dat$fcoi=='e'] <- 'he2e1A'
        catA[dat$alive & dat$ser==1 & dat$mcoi=='b' & dat$fcoi=='e'] <- 'hbeA'
        catA[dat$alive & dat$ser==1 & dat$mcoi=='e' & dat$fcoi=='b'] <- 'hebA'
        catA[dat$alive & dat$ser==1 & dat$mcoi=='b' & dat$fcoi=='p' & dat$fcoi.phase=='a'] <- 'hbpaA'
        catA[dat$alive & dat$ser==1 & dat$mcoi=='b' & dat$fcoi=='p' & !dat$fcoi.phase=='a'] <- 'hbpA'
        catA[dat$alive & dat$ser==1 & dat$mcoi=='p' & dat$fcoi=='b' & dat$mcoi.phase=='a'] <- 'hpbaA'
        catA[dat$alive & dat$ser==1 & dat$mcoi=='p' & dat$fcoi=='b' & !dat$mcoi.phase=='a'] <- 'hpbA'
        catA[dat$alive & dat$ser==1 & dat$mcoi=='e' & dat$fcoi=='p' & dat$fcoi.phase=='a'] <- 'hepaA'
        catA[dat$alive & dat$ser==1 & dat$mcoi=='e' & dat$fcoi=='p' & !dat$fcoi.phase=='a'] <- 'hepA'
        catA[dat$alive & dat$ser==1 & dat$mcoi=='p' & dat$fcoi=='e' & dat$mcoi.phase=='a'] <- 'hpeaA'
        catA[dat$alive & dat$ser==1 & dat$mcoi=='p' & dat$fcoi=='e' & !dat$mcoi.phase=='a'] <- 'hpeA'
        ## add to df
        dat$cat <- cat
        dat$catA <- catA
        return(dat)
      }


## compare fits to simulated data
compsim <- function(simdat,                # simulated data
                    parsout,               # outputted pars (including all PI's)
                    simpars,               # simulated parameters
                    dirnm,                 # where to save
                    browse = F)
  {
      if(browse) browser()
      ## index infections
      piGb1.sumI.lg <- simdat$cat %in% c('mb.', 'hbe', 'hbpa', 'hbp', 'hb1b2')
      piGe1.sumI.lg <- simdat$cat %in% c('me.', 'hepa', 'hep', 'he1e2')
      piG.b1sumI.lg <- simdat$cat %in% c('f.b', 'heb', 'hpba', 'hpb', 'hb2b1')
      piG.e1sumI.lg <- simdat$cat %in% c('f.e', 'hpea', 'hpe', 'he2e1')
      ind.infs <- sum(piGb1.sumI.lg) + sum(piGe1.sumI.lg) + sum(piG.b1sumI.lg) + sum(piG.e1sumI.lg)
      piGb1.sumI.r  <- sum(piGb1.sumI.lg) / ind.infs
      piGe1.sumI.r  <- sum(piGe1.sumI.lg) / ind.infs
      piG.b1sumI.r  <- sum(piG.b1sumI.lg) / ind.infs
      piG.e1sumI.r  <- sum(piG.e1sumI.lg) / ind.infs
###################################################################### 
      ## conditional on survival
      ## divide simulated data by routes of transmission for males
      simdatA <- simdat[simdat$alive,]
      pibUA.lg <- simdatA$catA %in% c("mb.A", "hbeA", "hbpaA", "hbpA", "hb1b2A", "hb2b1A")
      pieUA.lg <- simdatA$catA %in% c("me.A", "hepA", "hepaA", "hebA", "he1e2A", "he2e1A")
      pipUA.lg <- simdatA$catA %in% c("hpeaA", "hpeA", "hpbaA", "hpbA")
      pipUaA.lg <- simdatA$catA %in% c("hpeaA", "hpbaA")
      pipUcA.lg <- simdatA$catA %in% c("hpeA", "hpbA")    
      ## 
      inf.males <- sum(pibUA.lg) + sum(pieUA.lg) + sum(pipUA.lg)
      pibUA.r <- sum(pibUA.lg) / inf.males
      pieUA.r <- sum(pieUA.lg) / inf.males
      pipUA.r <- sum(pipUA.lg) / inf.males
      pipUaA.r <- sum(pipUaA.lg) / inf.males
      pipUcA.r <- sum(pipUcA.lg) / inf.males    
      ## divide simulated data by routes of transmission for females
      piUbA.lg <- simdatA$catA %in% c("f.bA", "hebA", "hpbaA", "hpbA", "hb1b2A", "hb2b1A")
      piUeA.lg <- simdatA$catA %in% c("f.eA", "hpeaA", "hpeA", "hbeA", "he1e2A", "he2e1A")
      piUpA.lg <- simdatA$catA %in% c("hepaA", "hbpaA", "hepA", "hbpA")
      piUpaA.lg <- simdatA$catA %in% c("hepaA", "hbpaA")
      piUpcA.lg <- simdatA$catA %in% c("hepA", "hbpA")    
      ## 
      inf.females <- sum(piUbA.lg) + sum(piUeA.lg) + sum(piUpA.lg)
      piUbA.r <- sum(piUbA.lg) / inf.females
      piUeA.r <- sum(piUeA.lg) / inf.females
      piUpA.r <- sum(piUpA.lg) / inf.females
      piUpaA.r <- sum(piUpaA.lg) / inf.females
      piUpcA.r <- sum(piUpcA.lg) / inf.females    
      ## Create a data frame that compares true routes of transmission
      ## for observed individuals to estimates
      bdis <- c(pibUA.r, pieUA.r, pipUA.r, pipUaA.r, pipUcA.r,
                piUbA.r, piUeA.r, piUpA.r, piUpaA.r, piUpcA.r,
                piGb1.sumI.r, piGe1.sumI.r, piG.b1sumI.r, piG.e1sumI.r)
      names(bdis) <- c("pibUA","pieUA","pipUA","pipUaA","pipUcA",
                       "piUbA","piUeA","piUpA","piUpaA","piUpcA",
                       "piGb1.sumI", "piGe1.sumI", "piG.b1sumI", "piG.e1sumI")
      ##
      hazs <- c("bmb","bfb","bme","bfe","bmp","bfp")
      show <- c(hazs, names(bdis))
    propbias <- data.frame(pars[show,], c(simpars[hazs],bdis))
    colnames(propbias) <- c("2.5%",    "50%",     "97.5%",   "trueval")
    propbias[hazs,] <- propbias[hazs,] * 12 * 100
    propbias$bias <- (propbias[,'50%'] - propbias[,'trueval']) / propbias[,'trueval']
    signif(propbias,2)
    ## not scaling for heterogeneity anymore because done inside simulation to keep geometric mean constant
    pdf(file.path(dirnm, "trans distr bias.pdf"), w = 8, h = 5)
                                        #    layout(matrix(c(1,2,7,10,3,4,8,11,5,6,9,12),4,3))
    par(mar = c(3,3.5,2,.5), oma = c(1,0,0,0)) # mfrow = c(3,2), 
    xlim <- c(0, max(propbias[,1:4]))
    plot(0,0, type = "n", bty = "n", ylab = "", yaxt = "n", ylim = c(.5, 7.5), 
         xlim = xlim)#,
                                        #main = paste(rownames(propbias)[jjj], ": bias =", signif(propbias$bias[jjj],2)))
    axis(2, at = 1:6, label = hazs, las = 2)
    for(jjj in 1:6)
      {
        points(propbias[jjj,'trueval'], jjj, pch = 19, col = "red", cex = 1.5)
        arrows(propbias[jjj,'2.5%'], jjj, propbias[jjj,'97.5%'], jjj, code = 3, angle = 90, length = .03)
        points(propbias[jjj,'50%'], jjj, pch = 19, cex = .8)
        if(jjj < 6) points(12*100*init.samp[,jjj], rep(jjj-.1, nrow(init.samp)), col = "blue", cex = .5, pch = 19)
        if(jjj==6) points(12*100*exp(init.samp[,6])*init.samp[,5], rep(jjj-.1, nrow(init.samp)), col = "blue", cex = .5, pch = 19)
      }
    legend("topright", c("true", "estimates (95% CIs)","chain initial conditions"), col = c("red","black","blue"), pch = 19, bty = "n")
    mtext('prevalence-standardized transmission hazard (per 100 person years)', side = 1, line = -.5, outer = T)
    dev.off()
    ## 
    pdf(file.path(dirnm, "pi distr bias.pdf"), w = 8, h = 3)
    par(mar = c(3,4.5,2,.5), mfrow = c(1,2))
    xlim <- c(0,1)
    for(ss in c(6,11))
      {
        plot(0,0, type = "n", bty = "n", ylab = "", yaxt = "n", ylim = c(.5, 5.5), 
             xlim = xlim)
        axis(2, at = 1:5, rownames(propbias)[(ss+1):(ss+5)], las = 2)
        for(jjj in 1:5)
          {
            points(propbias[jjj+ss,'trueval'], jjj, pch = 19, col = "red", cex = 1.5)
            arrows(propbias[jjj+ss,'2.5%'], jjj, propbias[jjj+ss,'97.5%'], jjj, code = 3, angle = 90, length = .03)
            points(propbias[jjj+ss,'50%'], jjj, pch = 19, cex = .8)
          }
      }
    legend("topright", c("true", "estimates (95% CIs)"), col = c("red","black"), pch = 19, bty = "n")
    dev.off()
    ## 
    pdf(file.path(dirnm, "pi index distr bias.pdf"), w = 5, h = 3)
    par(mar = c(3,5.5,2,.5))
    plot(0,0, type = "n", bty = "n", ylab = "", yaxt = "n", ylim = c(.5, 4.5), xlim = xlim)
    axis(2, at = 1:4, label = rownames(propbias)[17:20], las = 2)
    for(jjj in 17:20)
      {
        points(propbias[jjj,'trueval'], jjj-16, pch = 19, col = "red", cex = 1.5)
        arrows(propbias[jjj,'2.5%'], jjj-16, propbias[jjj,'97.5%'], jjj-16, code = 3, angle = 90, length = .03)
        points(propbias[jjj,'50%'], jjj-16, pch = 19, cex = .8)
      }
    legend("topright", c("true", "estimates (95% CIs)"), col = c("red","black"), pch = 19, bty = "n")
    dev.off()
  }
