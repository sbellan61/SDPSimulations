


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
