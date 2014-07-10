oldpcalc <- function(pars, dat, browse = F, compars = NULL, # compare lprob for true pars with other pars (fitted)
                  give.pis = F,              # return individual pi values for couples (for outside mcmc)
                  acute.sc = 7,                # acute phase infectiousness
                  sim = F,                   # just use this to simulate data given parameters? (outputs pser.a)
                  spar.log = F,                  # use sexual partner acquisition rate
                  heter = F,                     # individual heterogeneity; simulation only
                  ## using global variables for these
#                  partner.arv = F,               # ART coverage affects within-partnership transmission?
#                  low.coverage.arv = F,          # if T, only 50% of those on ART are not infectious
                  survive = T,               # account for survival in analysis
                  cond.sim = F,              # only simulate individuals that will live
                  lrho.sd = 1/2,             # sd on lrho prior
                  trace = T) # only do certain calculations when tracing parameters (i.e. for non-thinned versions)
  {
    if(low.coverage.arv)    # if assuming only 50% of those on ART are not-infectious
      {
        cov.scalar <- .5
      }else{
        cov.scalar <- 1
      }
    K <- nrow(dat)
    if(!sim)
      {
        hh.log <- dat$ser==1
        mm.log <- dat$ser==2
        ff.log <- dat$ser==3
        ss.log <- dat$ser==4
      }
    if(sum(pars[1:5]<0)>0)                   #if any parameters are <0 then the model must be rejected so we return logprob =-Inf
      {
        probs <- NA
        lprob <- -Inf
        pop.avs <- NA
        proj12 <- NA
        pser.a <- NA
        pser <- NA
        rrs <- NA
      }else{
        bmb <- as.numeric(pars["bmb"])
        bfb <- as.numeric(pars["bfb"])
        bme <- as.numeric(pars["bme"])
        bfe <- as.numeric(pars["bfe"])
        bmp <- as.numeric(pars["bmp"])
        rho <- exp(as.numeric(pars["lrho"])) # feeding in log(rho)
        bfp <- bmp * rho
        # L stands for *L*ast iteration
        s..L <- rep(1,K)                # s: concordant negative
        mb.a1L <- rep(0, K)               # m: male positive
        mb.a2L <- rep(0, K)               # m: male positive
        mb.L <- rep(0, K)               # m: male positive        
        me.a1L <- rep(0, K)
        me.a2L <- rep(0, K)
        me.L <- rep(0, K)                
        f.ba1L <- rep(0, K)               # f: female positive
        f.ba2L <- rep(0, K)               # f: female positive
        f.bL <- rep(0, K)               # f: female positive        
        f.ea1L <- rep(0, K)
        f.ea2L <- rep(0, K)
        f.eL <- rep(0, K)                
        hb1b2L <- rep(0, K)               # h: concordant positive
        hb2b1L <- rep(0, K)               # b1b2 is both inf before, but female 1st, & vice versa
        hbeL <- rep(0, K)
        hebL <- rep(0, K)               # 2nd & 3rd character give route of transmission for M & F, respectively
        hbpaL <- rep(0, K)              # acute
        hpbaL <- rep(0, K)
        hepaL <- rep(0, K)
        hpeaL <- rep(0, K)
        hbpL <- rep(0, K)               # chronic
        hpbL <- rep(0, K)        
        hepL <- rep(0, K)        
        hpeL <- rep(0, K)        
        he2e1L <- rep(0, K)
        he1e2L <- rep(0, K)
        ## initiate vectors to update based on *L*ast state
        s.. <- rep(1,K)                # s: concordant negative
        mb.a1 <- rep(0, K)               # m: male positive
        mb.a2 <- rep(0, K)               # m: male positive
        mb. <- rep(0, K)               # m: male positive        
        me.a1 <- rep(0, K)
        me.a2 <- rep(0, K)
        me. <- rep(0, K)                
        f.ba1 <- rep(0, K)               # f: female positive
        f.ba2 <- rep(0, K)               # f: female positive
        f.b <- rep(0, K)               # f: female positive        
        f.ea1 <- rep(0, K)
        f.ea2 <- rep(0, K)
        f.e <- rep(0, K)                
        hb1b2 <- rep(0, K)               # h: concordant positive
        hb2b1 <- rep(0, K)               # b1b2 is both inf before, but female 1st, & vice versa
        hbe <- rep(0, K)
        heb <- rep(0, K)               # 2nd & 3rd character give route of transmission for M & F, respectively
        hbpa <- rep(0, K)
        hpba <- rep(0, K)
        hepa <- rep(0, K)
        hpea <- rep(0, K)
        hbp <- rep(0, K)         
        hpb <- rep(0, K)        
        hep <- rep(0, K)        
        hpe <- rep(0, K)        
        he2e1 <- rep(0, K)
        he1e2 <- rep(0, K)
        # i.e., hbpa is a ++ couple in which the male was inf *b*efore
        # couple formation & the female by her *p*artner while he was *a*cutely infectious
        if(survive)
          {
            # A stands for *A*live, i.e. joint probability of serostatus
            # and both partners being alive at sampling
            mb.a1AL <- rep(0, K)               # m: male positive
            mb.a2AL <- rep(0, K)               # m: male positive
            mb.AL <- rep(0, K)               # m: male positive        
            me.a1AL <- rep(0, K)
            me.a2AL <- rep(0, K)
            me.AL <- rep(0, K)                
            f.ba1AL <- rep(0, K)               # f: female positive
            f.ba2AL <- rep(0, K)               # f: female positive
            f.bAL <- rep(0, K)               # f: female positive        
            f.ea1AL <- rep(0, K)
            f.ea2AL <- rep(0, K)
            f.eAL <- rep(0, K)                
            hb1b2AL <- rep(0, K)               # h: concordant positive
            hb2b1AL <- rep(0, K)               # b1b2 is both inf before, but female 1st, & vice versa
            hbeAL <- rep(0, K)
            hebAL <- rep(0, K)               # 2nd & 3rd character give route of transmission for M & F, respectively
            hbpaAL <- rep(0, K)
            hbpAL <- rep(0, K)         
            hpbaAL <- rep(0, K)
            hpbAL <- rep(0, K)        
            hepaAL <- rep(0, K)
            hepAL <- rep(0, K)        
            hpeaAL <- rep(0, K)
            hpeAL <- rep(0, K)        
            he2e1AL <- rep(0, K)
            he1e2AL <- rep(0, K)
            ## initiate vectors to update based on *L*ast state
            mb.a1A <- rep(0, K)               # m: male positive
            mb.a2A <- rep(0, K)               # m: male positive
            mb.A <- rep(0, K)               # m: male positive        
            me.a1A <- rep(0, K)
            me.a2A <- rep(0, K)
            me.A <- rep(0, K)                
            f.ba1A <- rep(0, K)               # f: female positive
            f.ba2A <- rep(0, K)               # f: female positive
            f.bA <- rep(0, K)               # f: female positive        
            f.ea1A <- rep(0, K)
            f.ea2A <- rep(0, K)
            f.eA <- rep(0, K)                
            hb1b2A <- rep(0, K)               # h: concordant positive
            hb2b1A <- rep(0, K)               # b1b2 is both inf before, but female 1st, & vice versa
            hbeA <- rep(0, K)
            hebA <- rep(0, K)               # 2nd & 3rd character give route of transmission for M & F, respectively
            hbpaA <- rep(0, K)
            hbpA <- rep(0, K)         
            hpbaA <- rep(0, K)
            hpbA <- rep(0, K)        
            hepaA <- rep(0, K)
            hepA <- rep(0, K)        
            hpeaA <- rep(0, K)
            hpeA <- rep(0, K)        
            he2e1A <- rep(0, K)
            he1e2A <- rep(0, K)                
          }
        if(browse) browser()
        for(tt in 1:max(dat$bd))
          {
            ## probabilities are non-zero only for times after started having sex and before couple formation
            m.sex <- dat$tmar-dat$bd+tt-1 >= dat$tms & dat$tmar-dat$bd+tt-1 < dat$tmar
            f.sex <- dat$tmar-dat$bd+tt-1 >= dat$tfs & dat$tmar-dat$bd+tt-1 < dat$tmar
            e.sex <- m.sex|f.sex           # either are active
            ## probability infected in month tt
            p.m.bef <- rep(0,K)
            p.f.bef <- rep(0,K)
            ## spar
            if(spar.log)
              {
                m.spar <- rep(1, K)
                f.spar <- rep(1, K)
              }else{
                m.spar <- rep(1, K)
                f.spar <- rep(1, K)
              }
            if(heter)
              {
                m.het <- exp( rnorm(K, mean = 0, sd = 1) )
                f.het <- exp( rnorm(K, mean = 0, sd = 1) )
              }else{
                m.het <- rep(1, K)
                f.het <- rep(1, K)
              }
            p.m.bef[m.sex] <- (1 - exp(-bmb * m.spar[m.sex] * m.het[m.sex] *epicf[cbind(dat$tmar[m.sex]-dat$bd[m.sex]+tt-1, dat$epic.ind[m.sex])]))
            p.f.bef[f.sex] <- (1 - exp(-bfb * f.spar[f.sex] * f.het[f.sex] *epicm[cbind(dat$tmar[f.sex]-dat$bd[f.sex]+tt-1, dat$epic.ind[f.sex])]))
            ## probability infected in month tt and alive at sampling
            if(survive)
              {
                p.m.bef.a <- rep(0,K)
                p.f.bef.a <- rep(0,K)
                ## csurv[time til interview, age in months in this month]
                p.m.bef.a[m.sex] <- p.m.bef[m.sex] * csurv[cbind(dat$mage[m.sex]-dat$cd[m.sex]-dat$bd[m.sex]+tt-1, dat$cd[m.sex]+dat$bd[m.sex]-tt+1)]
                p.f.bef.a[f.sex] <- p.f.bef[f.sex] * csurv[cbind(dat$fage[f.sex]-dat$cd[f.sex]-dat$bd[f.sex]+tt-1, dat$cd[f.sex]+dat$bd[f.sex]-tt+1)]
              }
            ## iterate probabilities based on previous values for only cases where it needs updating
            s..[e.sex]   <- s..L[e.sex]*(1-p.m.bef[e.sex])*(1-p.f.bef[e.sex])
            mb.a1[e.sex] <- s..L[e.sex]*p.m.bef[e.sex]*(1-p.f.bef[e.sex])
            mb.a2[e.sex] <- mb.a1L[e.sex]*(1 - p.f.bef[e.sex])
            mb.[e.sex]   <- mb.a2L[e.sex]*(1 - p.f.bef[e.sex]) + mb.L[e.sex]*(1 - p.f.bef[e.sex])
            f.ba1[e.sex] <- s..L[e.sex]*p.f.bef[e.sex]*(1-p.m.bef[e.sex])
            f.ba2[e.sex] <- f.ba1L[e.sex]*(1 - p.m.bef[e.sex])
            f.b[e.sex]   <- f.ba2L[e.sex]*(1 - p.m.bef[e.sex]) + f.bL[e.sex]*(1 - p.m.bef[e.sex])
            ## for individuals infected in the same month, assign
            ## the order of infection based on competing risks
            ## formula, but if the denominator is 0, replace both
            ## with 0 to avoid errors.
            p.mfirst <- p.m.bef[e.sex] / (p.m.bef[e.sex]+p.f.bef[e.sex])
            p.ffirst <- 1-p.mfirst
            p.mfirst[is.na(p.mfirst)] <- 0
            p.ffirst[is.na(p.ffirst)] <- 0                
            hb1b2[e.sex] <- hb1b2L[e.sex] + p.mfirst * s..L[e.sex]*p.m.bef[e.sex]*p.f.bef[e.sex] +
                                            (mb.a1L[e.sex] + mb.a2L[e.sex] + mb.L[e.sex]) * p.f.bef[e.sex]
            hb2b1[e.sex] <- hb2b1L[e.sex] + p.ffirst * s..L[e.sex]*p.m.bef[e.sex]*p.f.bef[e.sex] +
                                            (f.ba1L[e.sex] + f.ba2L[e.sex] + f.bL[e.sex]) * p.m.bef[e.sex]
            ## iterate joint probabilities with survival
            if(survive)
              {
                mb.a1A[e.sex] <- s..L[e.sex]*p.m.bef.a[e.sex]*(1-p.f.bef[e.sex])
                mb.a2A[e.sex] <- mb.a1AL[e.sex]*(1 - p.f.bef[e.sex])
                mb.A[e.sex]   <- mb.a2AL[e.sex]*(1 - p.f.bef[e.sex]) + mb.AL[e.sex]*(1 - p.f.bef[e.sex])
                f.ba1A[e.sex] <- s..L[e.sex]*p.f.bef.a[e.sex]*(1-p.m.bef[e.sex])
                f.ba2A[e.sex] <- f.ba1AL[e.sex]*(1 - p.m.bef[e.sex])
                f.bA[e.sex]   <- f.ba2AL[e.sex]*(1 - p.m.bef[e.sex]) + f.bAL[e.sex]*(1 - p.m.bef[e.sex])
                ## for individuals infected in the same month, assign
                ## the order of infection based on competing risks
                ## formula, but if the denominator is 0, replace both
                ## with 0 to avoid errors.
                p.mfirst.a <- p.m.bef.a[e.sex] / (p.m.bef.a[e.sex]+p.f.bef.a[e.sex])
                p.ffirst.a <- 1-p.mfirst.a
                p.mfirst.a[is.na(p.mfirst.a)] <- 0
                p.ffirst.a[is.na(p.ffirst.a)] <- 0                
                hb1b2A[e.sex] <- hb1b2AL[e.sex] + p.mfirst.a * s..L[e.sex]*p.m.bef.a[e.sex]*p.f.bef.a[e.sex] +
                                                  (mb.a1AL[e.sex] + mb.a2AL[e.sex] + mb.AL[e.sex]) * p.f.bef.a[e.sex]
                hb2b1A[e.sex] <- hb2b1AL[e.sex] + p.ffirst.a * s..L[e.sex]*p.m.bef.a[e.sex]*p.f.bef.a[e.sex] +
                                                  (f.ba1AL[e.sex] + f.ba2AL[e.sex] + f.bAL[e.sex]) * p.m.bef.a[e.sex]
                ## Update *A*live *L*ast states
                mb.a1AL[e.sex] <- mb.a1A[e.sex]
                mb.a2AL[e.sex] <- mb.a2A[e.sex]
                mb.AL[e.sex]   <- mb.A[e.sex]                
                f.ba1AL[e.sex] <- f.ba1A[e.sex]
                f.ba2AL[e.sex] <- f.ba2A[e.sex]
                f.bAL[e.sex]   <- f.bA[e.sex]                
                hb1b2AL[e.sex] <- hb1b2A[e.sex]
                hb2b1AL[e.sex] <- hb2b1A[e.sex]                
              }
            ## Update other *L*ast states
            s..L[e.sex]   <- s..[e.sex]
            mb.a1L[e.sex] <- mb.a1[e.sex]
            mb.a2L[e.sex] <- mb.a2[e.sex]
            mb.L[e.sex]   <- mb.[e.sex]            
            f.ba1L[e.sex] <- f.ba1[e.sex]
            f.ba2L[e.sex] <- f.ba2[e.sex]
            f.bL[e.sex]   <- f.b[e.sex]            
            hb1b2L[e.sex] <- hb1b2[e.sex]
            hb2b1L[e.sex] <- hb2b1[e.sex]                        
          }
            ## probability of being infected by partner (constant, used inside loop)
            p.m.part <- 1 - exp(-bmp)
            p.f.part <- 1 - exp(-bfp)
        ##  during acute
            p.m.part.ac <- 1 - exp(-acute.sc * bmp)
            p.f.part.ac <- 1 - exp(-acute.sc * bfp)

        allstates <- data.frame(s..,mb.a1,mb.a2,mb.,me.a1,me.a2,me.,f.ba1,f.ba2,f.b,f.ea1,f.ea2,f.e,
                                hb1b2,hb2b1,hbe,heb,
                                hbpa,hbp,
                                hpba,hpb,
                                hepa,hep,
                                hpea,hpe,
                                he1e2,he2e1,
                                mb.a1A,mb.a2A,mb.A,me.a1A,me.a2A,me.A,f.ba1A,f.ba2A,f.bA,f.ea1A,f.ea2A,f.eA,
                                hb1b2A,hb2b1A,hbeA,hebA,
                                hbpaA,hbpA,hpbaA,hpbA,hepaA,hepA,hpeaA,hpeA,
                                he1e2A,he2e1A)
        seros.pre.old <<- allstates

        ## Now loop through marriage
        for(tt in 1:max(dat$cd-1))
          {
            # browser()
            ## are partners formed in a couple?
            fmd <- dat$cd >= tt              
            ######################################################################
            ## everything below is automatically sum(fmd) length except p.m/f.part which are length 1
            ## probability infected extracouply in the ttc-th month of couple
            p.m.exc <- (1 - exp(-bme * m.spar[fmd] * m.het[fmd] *epicf[cbind(dat$tmar[fmd]+(tt-1), dat$epic.ind[fmd])]))
            p.f.exc <- (1 - exp(-bfe * f.spar[fmd] * f.het[fmd] *epicm[cbind(dat$tmar[fmd]+(tt-1), dat$epic.ind[fmd])]))
            ## adjust probability of being infected by partner by probability partner is infectious (i.e. not on ART)
            if(partner.arv)
              {
                p.m.part <- 1 - exp(-bmp * (1 - cov.scalar * art.prev[cbind(dat$tmar[fmd]+(tt-1), dat$epic.ind[fmd])]))
                p.f.part <- 1 - exp(-bfp * (1 - cov.scalar * art.prev[cbind(dat$tmar[fmd]+(tt-1), dat$epic.ind[fmd])]))
                ## acute
                p.m.part.ac <- 1 - exp(-acute.sc * bmp * (1 - cov.scalar * art.prev[cbind(dat$tmar[fmd]+(tt-1), dat$epic.ind[fmd])]))
                p.f.part.ac <- 1 - exp(-acute.sc * bfp * (1 - cov.scalar * art.prev[cbind(dat$tmar[fmd]+(tt-1), dat$epic.ind[fmd])]))
              }
            if(survive)
              {
                ## Survival probabilities
                s.p.m <- csurv[cbind(dat$mage[fmd]-dat$cd[fmd]+tt-1, dat$tint[fmd] - (dat$tmar[fmd] + tt - 1))]
                s.p.f <- csurv[cbind(dat$fage[fmd]-dat$cd[fmd]+tt-1, dat$tint[fmd] - (dat$tmar[fmd] + tt - 1))]
                ## Transmission probabilities from partner (jointly with survival)
                p.m.part.a <- p.m.part * s.p.m
                p.f.part.a <- p.f.part * s.p.f
                ## acute from partner (jointly with survival)
                p.m.part.a.ac <- p.m.part.ac * s.p.m
                p.f.part.a.ac <- p.f.part.ac * s.p.f
                ## extra-couple
                p.m.exc.a <- p.m.exc * s.p.m
                p.f.exc.a <- p.f.exc * s.p.f
              }
            ######################################################################
            ## iterate probabilities
            s..[fmd]   <- s..L[fmd]*(1-p.m.exc)*(1-p.f.exc)
            # the first month of marriage a male infected before marriage must be in stage a2 or later (note that if a male is infected in the month befor marriage, his partner still is exposed to 2 months of acute phase because all infection probabilitie depend on infection status in the last month)
            mb.a1[fmd] <- 0 
            mb.a2[fmd] <- mb.a1L[fmd]*(1-p.f.exc)*(1-p.f.part.ac)
            mb.[fmd]   <- mb.a2L[fmd]*(1-p.f.exc)*(1-p.f.part.ac) + mb.L[fmd]*(1-p.f.exc)*(1-p.f.part)            
            me.a1[fmd] <- s..L[fmd]*p.m.exc*(1-p.f.exc)
            me.a2[fmd] <- me.a1L[fmd]*(1-p.f.exc)*(1-p.f.part.ac) 
            me.[fmd]   <- me.a2L[fmd]*(1-p.f.exc)*(1-p.f.part.ac) + me.L[fmd]*(1-p.f.exc)*(1-p.f.part) 
            f.ba1[fmd] <- 0
            f.ba2[fmd] <- f.ba1L[fmd]*(1-p.m.exc)*(1-p.m.part.ac)
            f.b[fmd]   <- f.ba2L[fmd]*(1-p.m.exc)*(1-p.m.part.ac) + f.bL[fmd]*(1-p.m.exc)*(1-p.m.part)            
            f.ea1[fmd] <- s..L[fmd]*p.f.exc*(1-p.m.exc)
            f.ea2[fmd] <- f.ea1L[fmd]*(1-p.m.exc)*(1-p.m.part.ac)
            f.e[fmd]   <- f.ea2L[fmd]*(1-p.m.exc)*(1-p.m.part.ac) + f.eL[fmd]*(1-p.m.exc)*(1-p.m.part)
##          hb1b2[fmd] <- hb1b2L[fmd] # Doesn't change during couple duration
##          hb2b1[fmd] <- hb2b1L[fmd] # Doesn't change during couple duration                
            hbe[fmd]  <- hbeL[fmd]  + (mb.a1L[fmd] + mb.a2L[fmd])*(1-p.f.part.ac)*p.f.exc + mb.L[fmd]*(1-p.f.part)*p.f.exc
            heb[fmd]  <- hebL[fmd]  + (f.ba1L[fmd] + f.ba2L[fmd])*(1-p.m.part.ac)*p.m.exc + f.bL[fmd]*(1-p.m.part)*p.m.exc
            hbpa[fmd] <- hbpaL[fmd] + (mb.a1L[fmd] + mb.a2L[fmd])*p.f.part.ac
            hbp[fmd]  <- hbpL[fmd]  + mb.L[fmd]*p.f.part            
            hpba[fmd] <- hpbaL[fmd] + (f.ba1L[fmd] + f.ba2L[fmd])*p.m.part.ac
            hpb[fmd]  <- hpbL[fmd]  + f.bL[fmd]*p.m.part
            hepa[fmd] <- hepaL[fmd] + (me.a1L[fmd] + me.a2L[fmd])*p.f.part.ac
            hep[fmd]  <- hepL[fmd]  + me.L[fmd]*p.f.part            
            hpea[fmd] <- hpeaL[fmd] + (f.ea1L[fmd] + f.ea2L[fmd])*p.m.part.ac
            hpe[fmd]  <- hpeL[fmd]  + f.eL[fmd]*p.m.part            
            ## for individuals infected in the same month, assign
            ## the order of infection based on competing risks
            ## formula, but if the denominator is 0, replace both
            ## with 0 to avoid errors.
            p.mfirst <- p.m.exc / (p.m.exc+p.f.exc)
            p.ffirst <- 1-p.mfirst
            p.mfirst[is.na(p.mfirst)] <- 0
            p.ffirst[is.na(p.ffirst)] <- 0                
            he1e2[fmd] <- he1e2L[fmd] + p.mfirst * s..L[fmd]*p.m.exc*p.f.exc +
                                        (me.a1L[fmd] + me.a2L[fmd])*(1-p.f.part.ac)*p.f.exc +
                                        me.L[fmd]*(1-p.f.part)*p.f.exc
            he2e1[fmd] <- he2e1L[fmd] + p.ffirst * s..L[fmd]*p.m.exc*p.f.exc +
                                        (f.ea1L[fmd] + f.ea2L[fmd])*(1-p.m.part.ac)*p.m.exc +
                                        f.eL[fmd]*(1-p.m.part)*p.m.exc
            ######################################################################
            ## Iterate probabilities jointly with survival until survey.
            ## Note for probabilities of not being infected, we don't
            ## use the joint probability with being alive at sampling.
            if(survive)
              {
                mb.a1A[fmd] <- 0
                mb.a2A[fmd] <- mb.a1AL[fmd]*(1-p.f.exc)*(1-p.f.part.ac)
                mb.A[fmd]   <- mb.a2AL[fmd]*(1-p.f.exc)*(1-p.f.part.ac) + mb.AL[fmd]*(1-p.f.exc)*(1-p.f.part) 
                me.a1A[fmd] <- s..L[fmd]*p.m.exc.a*(1-p.f.exc)
                me.a2A[fmd] <- me.a1AL[fmd]*(1-p.f.exc)*(1-p.f.part.ac)
                me.A[fmd]   <- me.a2AL[fmd]*(1-p.f.exc)*(1-p.f.part.ac) + me.AL[fmd]*(1-p.f.exc)*(1-p.f.part)
                f.ba1A[fmd] <- 0
                f.ba2A[fmd] <- f.ba1AL[fmd]*(1-p.m.exc)*(1-p.m.part.ac)
                f.bA[fmd]   <- f.ba2AL[fmd]*(1-p.m.exc)*(1-p.m.part.ac) + f.bAL[fmd]*(1-p.m.exc)*(1-p.m.part)
                f.ea1A[fmd] <- s..L[fmd]*p.f.exc.a*(1-p.m.exc)
                f.ea2A[fmd] <- f.ea1AL[fmd]*(1-p.m.exc)*(1-p.m.part.ac)
                f.eA[fmd]   <- f.ea2AL[fmd]*(1-p.m.exc)*(1-p.m.part.ac) + f.eAL[fmd]*(1-p.m.exc)*(1-p.m.part)
##              hb1b2A[fmd] <- hb1b2AL[fmd] # Doesn't change during couple duration
##              hb2b1A[fmd] <- hb2b1AL[fmd] # Doesn't change during couple duration                
                hbeA[fmd]  <- hbeAL[fmd]  + (mb.a1AL[fmd] + mb.a2AL[fmd])*(1-p.f.part.ac)*p.f.exc.a + mb.AL[fmd]*(1-p.f.part)*p.f.exc.a
                hebA[fmd]  <- hebAL[fmd]  + (f.ba1AL[fmd] + f.ba2AL[fmd])*(1-p.m.part.ac)*p.m.exc.a + f.bAL[fmd]*(1-p.m.part)*p.m.exc.a
                hbpaA[fmd] <- hbpaAL[fmd] + (mb.a1AL[fmd] + mb.a2AL[fmd])*p.f.part.a.ac
                hbpA[fmd]  <- hbpAL[fmd]  + mb.AL[fmd]*p.f.part.a                
                hpbaA[fmd] <- hpbaAL[fmd] + (f.ba1AL[fmd] + f.ba2AL[fmd])*p.m.part.a.ac
                hpbA[fmd]  <- hpbAL[fmd]  + f.bAL[fmd]*p.m.part.a                
                hepaA[fmd] <- hepaAL[fmd] + (me.a1AL[fmd] + me.a2AL[fmd])*p.f.part.a.ac
                hepA[fmd]  <- hepAL[fmd]  + me.AL[fmd]*p.f.part.a                
                hpeaA[fmd] <- hpeaAL[fmd] + (f.ea1AL[fmd] + f.ea2AL[fmd])*p.m.part.a.ac
                hpeA[fmd]  <- hpeAL[fmd]  + f.eAL[fmd]*p.m.part.a                
                ## for individuals infected in the same month, assign
                ## the order of infection based on competing risks
                ## formula, but if the denominator is 0, replace both
                ## with 0 to avoid errors.
                p.mfirst.a <- p.m.exc.a / (p.m.exc.a+p.f.exc.a)
                p.ffirst.a <- 1-p.mfirst.a
                p.mfirst.a[is.na(p.mfirst.a)] <- 0
                p.ffirst.a[is.na(p.ffirst.a)] <- 0                
                he1e2A[fmd] <- he1e2AL[fmd] + p.mfirst.a * s..L[fmd]*p.m.exc.a*p.f.exc.a +
                                          (me.a1AL[fmd] + me.a2AL[fmd])*(1-p.f.part.ac)*p.f.exc.a +
                                          me.AL[fmd]*(1-p.f.part)*p.f.exc.a
                he2e1A[fmd] <- he2e1AL[fmd] + p.ffirst.a * s..L[fmd]*p.m.exc.a*p.f.exc.a +
                                          (f.ea1AL[fmd] + f.ea2AL[fmd])*(1-p.m.part.ac)*p.m.exc.a +
                                          f.eAL[fmd]*(1-p.m.part)*p.m.exc.a
                ## update *L*ast month states for *A*live states
                mb.a1AL[fmd] <- mb.a1A[fmd]
                mb.a2AL[fmd] <- mb.a2A[fmd]
                mb.AL[fmd]   <- mb.A[fmd]                
                me.a1AL[fmd] <- me.a1A[fmd]
                me.a2AL[fmd] <- me.a2A[fmd]
                me.AL[fmd]   <- me.A[fmd]                
                f.ba1AL[fmd] <- f.ba1A[fmd]
                f.ba2AL[fmd] <- f.ba2A[fmd]
                f.bAL[fmd]   <- f.bA[fmd]                
                f.ea1AL[fmd] <- f.ea1A[fmd]
                f.ea2AL[fmd] <- f.ea2A[fmd]
                f.eAL[fmd]   <- f.eA[fmd]                
##              hb1b2AL[fmd] <- hb1b2A[fmd]
##              hb2b1AL[fmd] <- hb2b1A[fmd]                
                hbeAL[fmd] <- hbeA[fmd]
                hebAL[fmd] <- hebA[fmd]
                hbpaAL[fmd] <- hbpaA[fmd]
                hbpAL[fmd] <- hbpA[fmd]                
                hpbaAL[fmd] <- hpbaA[fmd]
                hpbAL[fmd] <- hpbA[fmd]                
                hepaAL[fmd] <- hepaA[fmd]
                hepAL[fmd] <- hepA[fmd]                
                hpeaAL[fmd] <- hpeaA[fmd]
                hpeAL[fmd] <- hpeA[fmd]                              
                he1e2AL[fmd] <- he1e2A[fmd]
                he2e1AL[fmd] <- he2e1A[fmd]                              
              }
            ## update other *L*ast month states
            s..L[fmd]   <-  s..[fmd]
            mb.a1L[fmd] <-  mb.a1[fmd]
            mb.a2L[fmd] <-  mb.a2[fmd]
            mb.L[fmd]   <-  mb.[fmd]            
            me.a1L[fmd] <-  me.a1[fmd]
            me.a2L[fmd] <-  me.a2[fmd]
            me.L[fmd]   <-  me.[fmd]            
            f.ba1L[fmd] <-  f.ba1[fmd]
            f.ba2L[fmd] <-  f.ba2[fmd]
            f.bL[fmd]   <-  f.b[fmd]            
            f.ea1L[fmd] <-  f.ea1[fmd]
            f.ea2L[fmd] <-  f.ea2[fmd]
            f.eL[fmd]   <-  f.e[fmd]            
##          hb1b2L[fmd] <-  hb1b2[fmd]
##          hb2b1L[fmd] <-  hb2b1[fmd]            
            hbeL[fmd] <-  hbe[fmd]
            hebL[fmd] <-  heb[fmd]
            hbpaL[fmd] <-  hbpa[fmd]
            hbpL[fmd] <-  hbp[fmd]            
            hpbaL[fmd] <-  hpba[fmd]
            hpbL[fmd] <-  hpb[fmd]            
            hepaL[fmd] <-  hepa[fmd]
            hepL[fmd] <-  hep[fmd]            
            hpeaL[fmd] <-  hpea[fmd]
            hpeL[fmd] <-  hpe[fmd]                     
            he1e2L[fmd] <-  he1e2[fmd]
            he2e1L[fmd] <-  he2e1[fmd]

            allstates <- data.frame(s..,mb.a1,mb.a2,mb.,me.a1,me.a2,me.,f.ba1,f.ba2,f.b,f.ea1,f.ea2,f.e,
                                hb1b2,hb2b1,hbe,heb,
                                hbpa,hbp,
                                hpba,hpb,
                                hepa,hep,
                                hpea,hpe,
                                he1e2,he2e1,
                                mb.a1A,mb.a2A,mb.A,me.a1A,me.a2A,me.A,f.ba1A,f.ba2A,f.bA,f.ea1A,f.ea2A,f.eA,
                                hb1b2A,hb2b1A,hbeA,hebA,
                                hbpaA,hbpA,hpbaA,hpbA,hepaA,hepA,hpeaA,hpeA,
                                he1e2A,he2e1A)
            if(tt == 1) print(p.f.exc.a)
            if(tt==1) seros1.old <<- allstates
          }
######################################################################         
        allstates <- data.frame(s..,mb.a1,mb.a2,mb.,me.a1,me.a2,me.,f.ba1,f.ba2,f.b,f.ea1,f.ea2,f.e,
                                hb1b2,hb2b1,hbe,heb,
                                hbpa,hbp,
                                hpba,hpb,
                                hepa,hep,
                                hpea,hpe,
                                he1e2,he2e1,
                                mb.a1A,mb.a2A,mb.A,me.a1A,me.a2A,me.A,f.ba1A,f.ba2A,f.bA,f.ea1A,f.ea2A,f.eA,
                                hb1b2A,hb2b1A,hbeA,hebA,
                                hbpaA,hbpA,hpbaA,hpbA,hepaA,hepA,hpeaA,hpeA,
                                he1e2A,he2e1A)
        seros.post.old <<- allstates
        
        ss <- s..
        mm <- mb.a1 + me.a1 + mb.a2 + me.a2 + mb. + me.
        ff <- f.ba1 + f.ea1 + f.ba2 + f.ea2 + f.b + f.e
        hh <- hb1b2 + hb2b1 + hbe + heb + hbpa + hpba + hepa + hpea + hbp + hpb + hep + hpe + he1e2 + he2e1
        ## Calculate probability of data given parameters * priors of paramters
        if(survive)
          {
            mmA <- mb.a1A + me.a1A + mb.a2A + me.a2A + mb.A + me.A
            ffA <- f.ba1A + f.ea1A + f.ba2A + f.ea2A + f.bA + f.eA
            hhA <- hb1b2A + hb2b1A + hbeA + hebA + hbpaA + hpbaA + hepaA + hpeaA + hbpA + hpbA + hepA + hpeA + he1e2A + he2e1A
            pser.a <- cbind(hhA, mmA, ffA, ss)
            pser.a.old <<- pser.a
          }
        pser <- cbind(hh, mm, ff, ss)
        if(trace & sim) # calculate expected route of transmission breakdowns for couples with unknown (or simulated serostatus)
          {
######################################################################
            ## Route of transmission breakdowns for observed couples
            ## (conditional on survival)
######################################################################
            ## male breakdown amongst observed M+F- couples, partner *N*egative
            pibNA <- mean((mb.a1A + mb.a2A + mb.A) / mmA, na.rm = T) 
            pieNA <- mean((me.a1A + me.a2A + me.A) / mmA, na.rm = T) 
            ## female breakdown amongst observed M-F+ couples, partner *N*egative
            piNbA <- mean((f.ba1A + f.ba2A + f.bA) / ffA, na.rm = T)
            piNeA <- mean((f.ea1A + f.ea2A + f.eA) / ffA, na.rm = T)
            ## male breakdown amongst observed M+F+ couples,  partner *P*ositive
            pibPA <- mean((hb1b2A + hb2b1A + hbeA + hbpaA + hbpA) / hhA, na.rm = T) 
            piePA <- mean((hebA + hepaA + hepA + he1e2A + he2e1A) / hhA, na.rm = T) 
            pipPA <- mean((hpbaA + hpeaA + hpbA + hpeA) / hhA, na.rm = T)
            pipPaA <- mean((hpbaA + hpeaA) / hhA, na.rm = T) # males infected by their partner in M+F+ couples during her acute phase.
            pipPcA <- mean((hpbA + hpeA) / hhA, na.rm = T) # chronic phase
            ## female breakdown amongst observed M+F+ couples,  partner *P*ositive
            piPbA <- mean((hb1b2A + hb2b1A + hebA + hpbaA + hpbA) / hhA, na.rm = T) 
            piPeA <- mean((hbeA + hpeaA + hpeA + he1e2A + he2e1A) / hhA, na.rm = T) 
            piPpA <- mean((hbpaA + hepaA + hbpA + hepA) / hhA, na.rm = T)
            piPpaA <- mean((hbpaA + hepaA) / hhA, na.rm = T) # acute phase
            piPpcA <- mean((hbpA + hepA) / hhA, na.rm = T) # chronic
            ## male breakdown amongst infected males in any observed couples,  partner *U*nknown (bc could be either)
            pibUA <- mean((mb.a1A + mb.a2A + mb.A + hb1b2A + hb2b1A + hbeA + hbpaA + hbpA) / (mmA + hhA), na.rm = T)
            pieUA <- mean((me.a1A + me.a2A + me.A + hebA + hepaA + hepA + he1e2A + he2e1A) / (mmA + hhA), na.rm = T)
            pipUA <- mean((hpbaA + hpeaA + hpbA + hpeA) / (mmA + hhA), na.rm = T)
            pipUaA <- mean((hpbaA + hpeaA) / (mmA + hhA), na.rm = T) #acute
            pipUcA <- mean((hpbA + hpeA) / (mmA + hhA), na.rm = T) #chronic
            ## female breakdown amongst infected females in any observed couples,  partner *U*nknown (bc could be either)
            piUbA <- mean((f.ba1A + f.ba2A + f.bA + hb1b2A + hb2b1A + hebA + hpbaA + hpbA) / (ffA + hhA), na.rm = T)
            piUeA <- mean((f.ea1A + f.ea2A + f.eA + hbeA + hpeaA + hpeA + he1e2A + he2e1A) / (ffA + hhA), na.rm = T)
            piUpA <- mean((hbpaA + hepaA + hbpA + hepA) / (ffA + hhA), na.rm = T)
            piUpaA <- mean((hbpaA + hepaA) / (ffA + hhA), na.rm = T) #acute
            piUpcA <- mean((hbpA + hepA) / (ffA + hhA), na.rm = T) #chronic
######################################################################
            pop.avs <- data.frame(
                                  ## conditional on survival
                                  pibNA, pieNA, # b/e in +- given A
                                  piNbA, piNeA, # b/e in -+ given A
                                  pibPA, piePA, pipPA, pipPaA, pipPcA, # b/e/p in male in ++ given A
                                  piPbA, piPeA, piPpA, piPpaA, piPpcA, # b/e/p in female in ++ given A
                                  pibUA, pieUA, pipUA, pipUaA, pipUcA, # b/e/p in male in any given A
                                  piUbA, piUeA, piUpA, piUpaA, piUpcA) # b/e/p in female in any given A
          }
        if(trace & !sim)
          {
            ######################################################################
            ## Route of transmission breakdowns for observed couples
            ## (conditional on survival)
            ######################################################################
            ## male breakdown amongst observed M+F- couples, partner *N*egative
            pibNA <- sum((mb.a1A[mm.log] + mb.a2A[mm.log] + mb.A[mm.log]) / mmA[mm.log]) / sum(mm.log)
            pieNA <- sum((me.a1A[mm.log] + me.a2A[mm.log] + me.A[mm.log]) / mmA[mm.log]) / sum(mm.log)
            ## female breakdown amongst observed M-F+ couples, partner *N*egative
            piNbA <- sum((f.ba1A[ff.log] + f.ba2A[ff.log] + f.bA[ff.log]) / ffA[ff.log]) / sum(ff.log)
            piNeA <- sum((f.ea1A[ff.log] + f.ea2A[ff.log] + f.eA[ff.log]) / ffA[ff.log]) / sum(ff.log)
            ## male breakdown amongst observed M+F+ couples,  partner *P*ositive
            pibPA <- sum((hb1b2A[hh.log] + hb2b1A[hh.log] + hbeA[hh.log] + hbpaA[hh.log] + hbpA[hh.log]) / hhA[hh.log]) / sum(hh.log)
            piePA <- sum((hebA[hh.log] + hepaA[hh.log] + hepA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]) / sum(hh.log)
            pipPA <- sum((hpbaA[hh.log] + hpeaA[hh.log] + hpbA[hh.log] + hpeA[hh.log]) / hhA[hh.log]) / sum(hh.log)
            pipPaA <- sum((hpbaA[hh.log] + hpeaA[hh.log]) / hhA[hh.log]) / sum(hh.log) #acute
            pipPcA <- sum((hpbA[hh.log] + hpeA[hh.log]) / hhA[hh.log]) / sum(hh.log) #chronic
            ## female breakdown amongst observed M+F+ couples,  partner *P*ositive
            piPbA <- sum((hb1b2A[hh.log] + hb2b1A[hh.log] + hebA[hh.log] + hpbaA[hh.log] + hpbA[hh.log]) / hhA[hh.log]) / sum(hh.log)
            piPeA <- sum((hbeA[hh.log] + hpeaA[hh.log] + hpeA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]) / sum(hh.log)
            piPpA <- sum((hbpaA[hh.log] + hepaA[hh.log] + hbpA[hh.log] + hepA[hh.log]) / hhA[hh.log]) / sum(hh.log)
            piPpaA <- sum((hbpaA[hh.log] + hepaA[hh.log]) / hhA[hh.log]) / sum(hh.log) #acute
            piPpcA <- sum((hbpA[hh.log] + hepA[hh.log]) / hhA[hh.log]) / sum(hh.log) #chronic
            ## male breakdown amongst infected males in any observed couples,  partner *U*nknown (bc could be either)
            pibUA <- (pibNA*sum(mm.log) + pibPA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
            pieUA <- (pieNA*sum(mm.log) + piePA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
            pipUA <- (pipPA*sum(hh.log)) / (sum(mm.log) + sum(hh.log))
            pipUaA <- (pipPaA*sum(hh.log)) / (sum(mm.log) + sum(hh.log)) #acute, keep in mind this is only proportion of all M+ acutely infected by their *current* partner
            pipUcA <- (pipPcA*sum(hh.log)) / (sum(mm.log) + sum(hh.log)) #chronic, ditto above
            ## female breakdown amongst infected females in any observed couples,  partner *U*nknown (bc could be either)
            piUbA <- (piNbA*sum(ff.log) + piPbA*sum(hh.log)) / (sum(ff.log) + sum(hh.log))
            piUeA <- (piNeA*sum(ff.log) + piPeA*sum(hh.log)) / (sum(ff.log) + sum(hh.log))
            piUpA <- (piPpA*sum(hh.log)) / (sum(ff.log) + sum(hh.log))
            piUpaA <- (piPpaA*sum(hh.log)) / (sum(ff.log) + sum(hh.log)) #acute, see notes aboe
            piUpcA <- (piPpcA*sum(hh.log)) / (sum(ff.log) + sum(hh.log)) #chronic
            ######################################################################
            ## Give pieNA, piNeA, piePA, piPeA, for each couple
            if(give.pis)
              {
                ## probability infection was extracouple given ser
                piCe.A <- rep(NA, K)
                piC.eA <- rep(NA, K)
                piCe.A[mm.log] <- (me.a1A[mm.log] + me.a2A[mm.log] + me.A[mm.log]) / mmA[mm.log]
                piCe.A[hh.log] <- (hebA[hh.log] + hepaA[hh.log] + hepA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]
                piC.eA[ff.log] <- (f.ea1A[ff.log] + f.ea2A[ff.log] + f.eA[ff.log]) / ffA[ff.log]
                piC.eA[hh.log] <- (hbeA[hh.log] + hpeaA[hh.log] + hpeA[hh.log] + he1e2A[hh.log] + he2e1A[hh.log]) / hhA[hh.log]
                pis <- data.frame(piCe.A, piC.eA)
              }
            ######################################################################
            ## Route of transmission breakdowns for inferred
            ## pseudopopulation (unconditional on survival)
            ######################################################################
            ######################################################################
            ## Index infections, do with estimators summing over all
            ## infected couples and over all couples
            ## version 1 - all infected couples
            mb1. <- rowSums(allstates[!ss.log, c("mb.a1","mb.a2","mb.","hb1b2","hbe","hbpa","hbp")])
            me1. <- rowSums(allstates[!ss.log, c("me.a1","me.a2","me.","he1e2","hepa","hep")])
            ## female
            f.b1 <- rowSums(allstates[!ss.log, c("f.ba1","f.ba2","f.b","hb2b1","heb","hpba","hpb")])
            f.e1 <- rowSums(allstates[!ss.log, c("f.ea1","f.ea2","f.e","he2e1","hpea","hpe")])
            all.infA <- rowSums(allstates[!ss.log,names(allstates)[(grepl("m", names(allstates)) | grepl("f", names(allstates)) | grepl("h", names(allstates))) & grepl("A", names(allstates))]])
            ## number of inflated male before-couple index infections (in any couple)
            mb1.Infl <- sum( mb1. / all.infA)
            ## number of inflated male extra-couple index infections (in any couple)
            me1.Infl <- sum( me1. / all.infA)
            ## nufber of inflated male before-couple index infections (in any couple)
            f.b1Infl <- sum( f.b1 / all.infA)
            ## number of inflated male extra-couple index infections (in any couple)
            f.e1Infl <- sum( f.e1 / all.infA)
            ## number of inflated index infections (1 in each couple with an infection)
            IndInfl <- mb1.Infl + me1.Infl + f.b1Infl + f.e1Infl # note this includes h-classes..
            ######################################################################
            ## Proportion of index infections pooling across gender
            piGb1.sumI <- mb1.Infl / IndInfl
            piGe1.sumI <- me1.Infl / IndInfl
            piG.b1sumI <- f.b1Infl / IndInfl
            piG.e1sumI <- f.e1Infl / IndInfl
            ## ## version 2 - sum over all couples
            mb1. <- rowSums(allstates[, c("mb.a1","mb.a2","mb.","hb1b2","hbe","hbpa","hbp")])
            me1. <- rowSums(allstates[, c("me.a1","me.a2","me.","he1e2","hepa","hep")])
            ## female
            f.b1 <- rowSums(allstates[, c("f.ba1","f.ba2","f.b","hb2b1","heb","hpba","hpb")])
            f.e1 <- rowSums(allstates[, c("f.ea1","f.ea2","f.e","he2e1","hpea","hpe")])
            all.infA <- rowSums(allstates[,c("s..", names(allstates)[(grepl("m", names(allstates)) | grepl("f", names(allstates)) | grepl("h", names(allstates))) & grepl("A", names(allstates))])])
            ## number of inflated male before-couple index infections (in any couple)
            mb1.Infl <- sum( mb1. / all.infA)
            ## number of inflated male extra-couple index infections (in any couple)
            me1.Infl <- sum( me1. / all.infA)
            ## nufber of inflated male before-couple index infections (in any couple)
            f.b1Infl <- sum( f.b1 / all.infA)
            ## number of inflated male extra-couple index infections (in any couple)
            f.e1Infl <- sum( f.e1 / all.infA)
            ## number of inflated index infections (1 in each couple with an infection)
            IndInfl <- mb1.Infl + me1.Infl + f.b1Infl + f.e1Infl
######################################################################
            ## Proportion of index infections pooling across gender
            piGb1.sumIS <- mb1.Infl / IndInfl
            piGe1.sumIS <- me1.Infl / IndInfl
            piG.b1sumIS <- f.b1Infl / IndInfl
            piG.e1sumIS <- f.e1Infl / IndInfl
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
            ######################################################################
            ## Project incidence forward 12 months for each of the 3
            ## couple types (ss, mm, ff) for each country in the data
            ## set (because they have different population prevalences
            num.country <- length(unique(dat$epic.ind))
            cc.inds <- unique(dat$epic.ind)
            ## concordant negative
            ss12.ssL <- rep(1, num.country)
            mm12a1.ssL <- rep(0, num.country)
            mm12a2.ssL <- rep(0, num.country)
            mm12.ssL <- rep(0, num.country)            
            ff12a1.ssL <- rep(0, num.country)
            ff12a2.ssL <- rep(0, num.country)
            ff12.ssL <- rep(0, num.country)            
            hh12.ssL <- rep(0, num.country)            
            ## male positive discordant
            mm12a1.mmL <- rep(1, num.country)
            mm12a2.mmL <- rep(1, num.country)
            mm12.mmL <- rep(1, num.country)            
            hh12.mmL <- rep(0, num.country)            
            ## female positive discordant
            ff12a1.ffL <- rep(1, num.country)
            ff12a2.ffL <- rep(1, num.country)
            ff12.ffL <- rep(1, num.country)            
            hh12.ffL <- rep(0, num.country)
            ## initialize pis
            pi.m.part12.ss.ac <- 0      #acute
            pi.f.part12.ss.ac <- 0      #acute
            pi.m.part12.ss <- 0
            pi.f.part12.ss <- 0
            pi.m.exc12.ss <- 0
            pi.f.exc12.ss <- 0
            pi.f.part12.mm.ac <- 0      #acute
            pi.f.part12.mm <- 0            
            pi.f.exc12.mm <- 0            
            pi.m.part12.ff.ac <- 0      #acute
            pi.m.part12.ff <- 0            
            pi.m.exc12.ff <- 0            
            for(tt in 1:12)
              {
                if(partner.arv)         # put partner's on ART with some probability assigned related to ART coverage
                  {
                    p.m.part <- 1 - exp(-bmp * (1 - cov.scalar * art.prev[1332+tt-1, cc.inds]))
                    p.f.part <- 1 - exp(-bfp * (1 - cov.scalar * art.prev[1332+tt-1, cc.inds]))
                    p.m.part.ac <- 1 - exp(-acute.sc * bmp * (1 - cov.scalar * art.prev[1332+tt-1, cc.inds])) # acute
                    p.f.part.ac <- 1 - exp(-acute.sc * bfp * (1 - cov.scalar * art.prev[1332+tt-1, cc.inds])) # acute
                  }
                ######################################################################
                ## Transmission probabilities
                ## probability infected extracouply in various months of 2011
                p.m.exc <- 1 - exp(-bme*epicf[1332+tt-1, cc.inds])
                p.f.exc <- 1 - exp(-bfe*epicm[1332+tt-1, cc.inds])
                ## concordant negative couples
                ss12.ss   <- ss12.ssL*(1-p.m.exc)*(1-p.f.exc)
                mm12a1.ss <- ss12.ssL*p.m.exc*(1-p.f.exc)
                mm12a2.ss <- mm12a1.ssL*(1-p.f.exc)*(1-p.f.part.ac)
                mm12.ss   <- mm12a2.ssL*(1-p.f.exc)*(1-p.f.part.ac) + mm12.ssL*(1-p.f.exc)*(1-p.f.part)
                ff12a1.ss <- ss12.ssL*p.f.exc*(1-p.m.exc)
                ff12a2.ss <- ff12a1.ssL*(1-p.m.exc)*(1-p.m.part.ac)
                ff12.ss   <- ff12a2.ssL*(1-p.m.exc)*(1-p.m.part.ac) + ff12.ssL*(1-p.m.exc)*(1-p.m.part)
                hh12.ss <- hh12.ssL + ss12.ssL* p.m.exc*p.f.exc +
                    (mm12a1.ssL + mm12a2.ssL)*(p.f.part.ac + (1-p.f.part.ac)*p.f.exc) + mm12.ssL*(p.f.part + (1-p.f.part)*p.f.exc) +
                    (ff12a1.ssL + ff12a2.ssL)*(p.m.part.ac + (1-p.m.part.ac)*p.m.exc) + ff12.ssL*(p.m.part + (1-p.m.part)*p.m.exc)
                pi.m.part12.ss <- pi.m.part12.ss + (ff12a1.ssL + ff12a2.ssL)*p.m.part.ac + ff12.ssL*p.m.part
                pi.f.part12.ss <- pi.f.part12.ss + (mm12a1.ssL + mm12a2.ssL)*p.f.part.ac + mm12.ssL*p.f.part        
                pi.m.part12.ss.ac <- pi.m.part12.ss.ac + (ff12a1.ssL + ff12a2.ssL)*p.m.part.ac  # acute
                pi.f.part12.ss.ac <- pi.f.part12.ss.ac + (mm12a1.ssL + mm12a2.ssL)*p.f.part.ac  # acute
                pi.m.exc12.ss <- pi.m.exc12.ss + (ss12.ssL + (ff12a1.ssL + ff12a2.ssL)*(1-p.m.part.ac) + ff12.ssL*(1-p.m.part))*p.m.exc
                pi.f.exc12.ss <- pi.f.exc12.ss + (ss12.ssL + (mm12a1.ssL + mm12a2.ssL)*(1-p.f.part.ac) + mm12.ssL*(1-p.f.part))*p.f.exc
                ## male positive couples & female seroconversion,
                ##  assume there are no prevalent serodiscordant couples with individuals in acute stage at DHS ****
                mm12a1.mm <- 0          # all M+ couples  are already M+ at start of 12 month projection
                mm12a2.mm <- 0
                mm12.mm   <- mm12.mmL*(1-p.f.exc)*(1-p.f.part)                
                hh12.mm   <- hh12.mmL + mm12.mmL*(p.f.part + (1-p.f.part)*p.f.exc) 
                pi.f.part12.mm <- pi.f.part12.mm + mm12.mmL*p.f.part        
                pi.f.exc12.mm <- pi.f.exc12.mm + mm12.mmL*(1-p.f.part)*p.f.exc
                ## female positive couples & male seroconversion                  
                ff12a1.ff <- 0
                ff12a2.ff <- 0
                ff12.ff   <- ff12.ffL*(1-p.m.exc)*(1-p.m.part)                
                hh12.ff   <- hh12.ffL + ff12.ffL*(p.m.part + (1-p.m.part)*p.m.exc)
                pi.m.part12.ff <- pi.m.part12.ff + ff12.ffL*p.m.part
                pi.m.exc12.ff <- pi.m.exc12.ff + ff12.ffL*(1-p.m.part)*p.m.exc
                ss12.ssL   <- ss12.ss
                mm12a1.ssL <- mm12a1.ss
                mm12a2.ssL <- mm12a2.ss
                mm12.ssL   <- mm12.ss                
                ff12a1.ssL <- ff12a1.ss
                ff12a2.ssL <- ff12a2.ss
                ff12.ssL   <- ff12.ss                
                hh12.ssL   <- hh12.ss
                ## male positive discordant
                mm12a1.mmL <- mm12a1.mm
                mm12a2.mmL <- mm12a2.mm
                mm12.mmL   <- mm12.mm
                hh12.mmL   <- hh12.mm
                ## female positive discordant
                ff12a1.ffL <- ff12a1.ff
                ff12a2.ffL <- ff12a2.ff
                ff12.ffL   <- ff12.ff                
                hh12.ffL   <- hh12.ff
              }
            n.m.part.dc <- 0
            n.m.part.cc.ac <- 0
            n.m.part.cc <- 0            
            n.f.part.dc <- 0
            n.f.part.cc.ac <- 0
            n.f.part.cc <- 0            
            n.m.exc.dc <- 0
            n.m.exc.cc <- 0
            n.f.exc.dc <- 0
            n.f.exc.cc <- 0
            ## add all the different countries incidence by scaling by serotype
            for(cc in 1:num.country)
              {
                n.m.part.dc <- n.m.part.dc + pi.m.part12.ff[cc]*sum(ff.log & dat$epic.ind == cc.inds[cc])
                n.f.part.dc <- n.f.part.dc + pi.f.part12.mm[cc]*sum(mm.log & dat$epic.ind == cc.inds[cc])                
                n.m.part.cc <- n.m.part.cc + pi.m.part12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])
                n.f.part.cc <- n.f.part.cc + pi.f.part12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])                
                n.m.part.cc.ac <- n.m.part.cc.ac + pi.m.part12.ss.ac[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])
                n.f.part.cc.ac <- n.f.part.cc.ac + pi.f.part12.ss.ac[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])                
                n.m.exc.dc <- n.m.exc.dc + pi.m.exc12.ff[cc]*sum(ff.log & dat$epic.ind == cc.inds[cc])
                n.f.exc.dc <- n.f.exc.dc + pi.f.exc12.mm[cc]*sum(mm.log & dat$epic.ind == cc.inds[cc])                
                n.m.exc.cc <- n.m.exc.cc + pi.m.exc12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])
                n.f.exc.cc <- n.f.exc.cc + pi.f.exc12.ss[cc]*sum(ss.log & dat$epic.ind == cc.inds[cc])                
              }
            n.m.dc <- n.m.part.dc + n.m.exc.dc
            n.f.dc <- n.f.part.dc + n.f.exc.dc
            n.m.part.tot <- n.m.part.dc + n.m.part.cc
            n.f.part.tot <- n.f.part.dc + n.f.part.cc            
            n.m.exc.tot <- n.m.exc.dc + n.m.exc.cc
            n.f.exc.tot <- n.f.exc.dc + n.f.exc.cc
            proj12 <- data.frame(n.m.part.dc, n.f.part.dc, # incidence per 1000
                                 n.m.part.cc, n.f.part.cc,
                                 n.m.part.cc.ac, n.f.part.cc.ac,
                                 n.m.exc.dc, n.f.exc.dc,
                                 n.m.exc.cc, n.f.exc.cc,
                                 n.m.part.tot, n.f.part.tot,
                                 n.m.exc.tot, n.f.exc.tot) / sum(!hh.log) * 1000
            prop.exc.m <- n.m.exc.tot / (n.m.exc.tot + n.m.part.tot)
            prop.exc.f <- n.f.exc.tot / (n.f.exc.tot + n.f.part.tot)
            prop.exc.m.dc <- n.m.exc.dc / n.m.dc
            prop.exc.f.dc <- n.f.exc.dc / n.f.dc
            proj12 <- data.frame(proj12, prop.exc.m, prop.exc.f, prop.exc.m.dc, prop.exc.f.dc) 
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
            rr.mf.bef.cont <- rr.mf.bef * rho
            rr.mf.exc.cont <- rr.mf.exc * rho
            rrs <- data.frame(rr.m.eb = rr.m.eb, rr.f.eb = rr.f.eb,
                              rr.m.pe = rr.m.pe, rr.f.pe = rr.f.pe,
                              rr.m.pb = rr.m.pb, rr.f.pb = rr.f.pb,
                              rr.mf.bef = rr.mf.bef, rr.mf.exc = rr.mf.exc,
                              rr.mf.bef.cont = rr.mf.bef.cont, rr.mf.exc.cont = rr.mf.exc.cont)
          }
        if(sim) # if simulating data
          {
            probs <- NA
            lprob <- NA
            ## create couple state probability *A*live & *D*ead
            sim.probs <- data.frame(s..A = s..,
                                    mb.a1A,   mb.a1D   =  mb.a1 - mb.a1A,
                                    me.a1A,   me.a1D   =  me.a1 - me.a1A,
                                    f.ba1A,   f.ba1D   =  f.ba1 - f.ba1A,
                                    f.ea1A,   f.ea1D   =  f.ea1 - f.ea1A,
                                    mb.a2A,   mb.a2D   =  mb.a2 - mb.a2A,
                                    me.a2A,   me.a2D   =  me.a2 - me.a2A,
                                    f.ba2A,   f.ba2D   =  f.ba2 - f.ba2A,
                                    f.ea2A,   f.ea2D   =  f.ea2 - f.ea2A,
                                    mb.A,   mb.D   =  mb.- mb.A,
                                    me.A,   me.D   =  me. - me.A,
                                    f.bA,   f.bD   =  f.b - f.bA,
                                    f.eA,   f.eD   =  f.e - f.eA,
                                    hb1b2A, hb1b2D = hb1b2 - hb1b2A, # note some of the dead cases were infected by dead partners, so can only use the index cases in any h couple
                                    hb2b1A, hb2b1D = hb2b1 - hb2b1A,
                                    hbeA,   hbeD   =  hbe - hbeA,
                                    hebA,   hebD   =  heb - hebA,
                                    hepaA,   hepaD   =  hepa - hepaA,
                                    hpeaA,   hpeaD   =  hpea - hpeaA,
                                    hbpaA,   hbpaD   =  hbpa - hbpaA,
                                    hpbaA,   hpbaD   =  hpba - hpbaA,
                                    hepA,   hepD   =  hep - hepA,
                                    hpeA,   hpeD   =  hpe - hpeA,
                                    hbpA,   hbpD   =  hbp - hbpA,
                                    hpbA,   hpbD   =  hpb - hpbA,
                                    he1e2A, he1e2D = he1e2 - he1e2A,
                                    he2e1A, he2e1D = he2e1 - he2e1A)
            for(ii in 1:nrow(dat))
              {
                        dat$cat[ii] <- which(rmultinom(1, 1, sim.probs[ii,])==1)
              }
            dat$cat.nm <- names(sim.probs)[dat$cat]
            dat$cat.nm <- factor(dat$cat.nm, levels = names(sim.probs))
            K <- nrow(dat)
            if(!survive) pser.a <- NA
          }else{ ## if not simulating data calculate likelihood p(data|pars)
            if(survive)
              {                         # must NORMALIZE probabilities to 1 for likelihood!
                probs <- pser.a[cbind(1:K,dat$ser)] / rowSums(pser.a) # accounting for survival
              }else{
                probs <- pser[cbind(1:K,dat$ser)] # if not accounting for survival
                pser.a <- NA
              }
            if(sum(probs==0)==0) # if non of the serotatuses occur with 0 probability in the current model
              {
                lprob <- sum(log(probs)) + dnorm(log(rho), log(trans.ratio), lrho.sd, log = T)
                if(length(compars)>0) clprob <- sum(log(cprobs)) + dnorm(as.numeric(compars["lrho"]),
                                                                         log(trans.ratio), lrho.sd, log = T)
              }else{ # if some of the serostatuses are 0, then the current parameters have 0 probability
                lprob <- -Inf
              }
          }
      }
    if(length(compars)==0)
      {
        clprob <- NA
        cprobs <- NA
      }
    if(sim)
      {
        if(trace)
          {
            if(give.pis)
              {
                return(list(lprob = lprob,pop.avs = pop.avs, proj12 = proj12, sim.probs, allstates = allstates,
                            pser.a = pser.a, pser = pser, dat = dat, clprob = clprob, probs = probs, cprobs = cprobs, m.het = m.het, f.het = f.het))
              }else{ 
                return(list(lprob = lprob,pop.avs = pop.avs, rrs = rrs, proj12=proj12, sim.probs, allstates = allstates,
                            pser.a = pser.a, pser = pser, dat = dat, clprob = clprob, probs = probs, cprobs = cprobs, m.het = m.het, f.het = f.het))
              }
          }else{
            return(list(lprob = lprob, pser.a = pser.a, pser = pser, dat = dat, sim.probs, allstates = allstates,
                        clprob = clprob, probs = probs, cprobs = cprobs, m.het = m.het, f.het = f.het))
          }
      }else{                            # if not simulating
        if(trace)
          {
            if(give.pis)
              {
                return(list(lprob = lprob,pop.avs = pop.avs, rrs = rrs,  proj12=proj12, pis = pis, allstates = allstates,
                            pser.a = pser.a, pser = pser, probs = probs))
              }else{ 
                return(list(lprob = lprob,pop.avs = pop.avs, rrs = rrs, proj12=proj12,
                            pser.a = pser.a, pser = pser, probs = probs))
              }
          }else{
            return(list(lprob = lprob, pser.a = pser.a, pser = pser))
          }
      }
  }

