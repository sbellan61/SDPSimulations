library(data.table); library(compiler)

pre.coupleDT <- function(serostates, sexually.active) {
    serostates[sexually.active , `:=`(
        s..   = s.. * (1-p.m.bef) * (1-p.f.bef),
        mb.a1 = s.. * p.m.bef * (1-p.f.bef),
        mb.a2 = mb.a1 * (1 - p.f.bef),
        mb.   = mb.a2 * (1 - p.f.bef) + mb. * (1 - p.f.bef),
        f.ba1 = s.. * p.f.bef * (1-p.m.bef),
        f.ba2 = f.ba1 * (1 - p.m.bef),
        f.b   = f.ba2 * (1 - p.m.bef) + f.b * (1 - p.m.bef),
        hb1b2 = hb1b2 + .5  *  s.. * p.m.bef * p.f.bef + (mb.a1 + mb.a2 + mb.)  *  p.f.bef,
        hb2b1 = hb2b1 + .5  *  s.. * p.m.bef * p.f.bef + (f.ba1 + f.ba2 + f.b)  *  p.m.bef)
           ]
    return(serostates)
}



pre.coupleMat <- function(serostates, sexually.active) {
    temp <- serostates[sexually.active,]
    serostates[sexually.active,'s..']   = temp[,'s..'] * (1-p.m.bef) * (1-p.f.bef)
    serostates[sexually.active,'mb.a1'] = temp[,'s..'] * p.m.bef * (1-p.f.bef)
    serostates[sexually.active,'mb.a2'] = temp[,'mb.a1'] * (1 - p.f.bef)
    serostates[sexually.active,'mb.'] = temp[,'mb.a2'] * (1 - p.f.bef) + temp[,'mb.'] * (1 - p.f.bef)
    serostates[sexually.active,'f.ba1'] = temp[,'s..'] * p.f.bef * (1-p.m.bef)
    serostates[sexually.active,'f.ba2'] = temp[,'f.ba1'] * (1 - p.m.bef)
    serostates[sexually.active,'f.b'] = temp[,'f.ba2'] * (1 - p.m.bef) + temp[,'f.b'] * (1 - p.m.bef)
    serostates[sexually.active,'hb1b2'] = temp[,'hb1b2'] + .5  *  temp[,'s..'] * p.m.bef * p.f.bef + (temp[,'mb.a1'] + temp[,'mb.a2'] + temp[,'mb.'])  *  p.f.bef
    serostates[sexually.active,'hb2b1'] = temp[,'hb2b1'] + .5  *  temp[,'s..'] * p.m.bef * p.f.bef + (temp[,'f.ba1'] + temp[,'f.ba2'] + temp[,'f.b'])  *  p.m.bef
    return(serostates)
}

for(ii in 1:length(state.var.nms)) assign(paste0(state.var.nms[ii],'Ind'), ii)
pre.coupleMatInd <- function(serostates, sexually.active) {
    temp <- serostates[sexually.active,]
    serostates[sexually.active,s..Ind]   = temp[,s..Ind] * (1-p.m.bef) * (1-p.f.bef)
    serostates[sexually.active,mb.a1Ind] = temp[,s..Ind] * p.m.bef * (1-p.f.bef)
    serostates[sexually.active,mb.a2Ind] = temp[,mb.a1Ind] * (1 - p.f.bef)
    serostates[sexually.active,mb.Ind] = temp[,mb.a2Ind] * (1 - p.f.bef) + temp[,mb.Ind] * (1 - p.f.bef)
    serostates[sexually.active,f.ba1Ind] = temp[,s..Ind] * p.f.bef * (1-p.m.bef)
    serostates[sexually.active,f.ba2Ind] = temp[,f.ba1Ind] * (1 - p.m.bef)
    serostates[sexually.active,f.bInd] = temp[,f.ba2Ind] * (1 - p.m.bef) + temp[,f.bInd] * (1 - p.m.bef)
    serostates[sexually.active,hb1b2Ind] = temp[,hb1b2Ind] + .5  *  temp[,s..Ind] * p.m.bef * p.f.bef + (temp[,mb.a1Ind] + temp[,mb.a2Ind] + temp[,mb.Ind])  *  p.f.bef
    serostates[sexually.active,hb2b1Ind] = temp[,hb2b1Ind] + .5  *  temp[,s..Ind] * p.m.bef * p.f.bef + (temp[,f.ba1Ind] + temp[,f.ba2Ind] + temp[,f.bInd])  *  p.m.bef
    return(serostates)
}

pre.coupleMatOld <- function(serostatesList, e.sex) {
    serostatesList <- with(serostatesList,
                           {
                               for(ii in 1:length(state.var.nms)) assign(paste0(state.var.nms[ii],'L'), get(state.var.nms[ii]))
                               s..[e.sex]   <- s..L[e.sex]*(1-p.m.bef)*(1-p.f.bef)
                               mb.a1[e.sex] <- s..L[e.sex]*p.m.bef*(1-p.f.bef)
                               mb.a2[e.sex] <- mb.a1L[e.sex]*(1 - p.f.bef)
                               mb.[e.sex]   <- mb.a2L[e.sex]*(1 - p.f.bef) + mb.L[e.sex]*(1 - p.f.bef)
                               f.ba1[e.sex] <- s..L[e.sex]*p.f.bef*(1-p.m.bef)
                               f.ba2[e.sex] <- f.ba1L[e.sex]*(1 - p.m.bef)
                               f.b[e.sex]   <- f.ba2L[e.sex]*(1 - p.m.bef) + f.bL[e.sex]*(1 - p.m.bef)
                               ## for individuals infected in the same month, assign
                               ## the order of infection based on competing risks
                               ## formula, but if the denominator is 0, replace both
                               ## with 0 to avoid errors.
                               hb1b2[e.sex] <- hb1b2L[e.sex] + .5 * s..L[e.sex]*p.m.bef*p.f.bef +
                                   (mb.a1L[e.sex] + mb.a2L[e.sex] + mb.L[e.sex]) * p.f.bef
                               hb2b1[e.sex] <- hb2b1L[e.sex] + .5 * s..L[e.sex]*p.m.bef*p.f.bef +
                                   (f.ba1L[e.sex] + f.ba2L[e.sex] + f.bL[e.sex]) * p.m.bef
                               return(serostatesList)
                           }
                           )
    return(serostates)
}

state.var.nms <- c('s..', 'mb.a1', 'mb.a2', 'mb.', 'f.ba1', 'f.ba2', 'f.b', 'hb1b2', 'hb2b1')
n <- 10^4
k <- 9
serostates <- matrix(0,n,k)
serostates <- as.data.table(serostates)
setnames(serostates, 1:k, state.var.nms)
serostates[, `:=`(s.. = 1)]
serostates
serostatesMat <- as.matrix(serostates)
serostatesMatInd <- serostatesMat
zros <- rep(0, n)
for(ii in 2:length(state.var.nms)) assign(state.var.nms[ii], zros)
s.. <- rep(1,n)
serostatesList <- list(s.., mb.a1, mb.a2, mb., f.ba1, f.ba2, f.b, hb1b2, hb2b1)
names(serostatesList) <- state.var.nms
length(serostatesList)
head(serostatesList[['s..']],5)
head(serostatesList[['hb1b2']],5)
p.m.bef <- .012
p.f.bef <- .07
its <- 12*40
sexually.activeMat <- matrix(rbinom(n*its, 1,.5),n,its)==1

tDT <- system.time(
    for(ii in 1:its) {
        serostates <- pre.coupleDT(serostates, sexually.activeMat[,ii])
    }
    ) ## about 2.25 seconds

tMat <- system.time(
    for(ii in 1:its) {
        serostatesMat <- pre.coupleMat(serostatesMat, sexually.activeMat[,ii])
    }
    ) ## about 6 seconds

tMatInd <- system.time(
    for(ii in 1:its) {
        serostatesMatInd <- pre.coupleMatInd(serostatesMatInd, sexually.activeMat[,ii])
    }
    ) ## about 6 seconds
tList <- system.time(
    for(ii in 1:its) {
        serostatesList <- pre.coupleMatOld(serostatesList, sexually.activeMat[,ii])
    }
    ) ## about 9 seconds
## tglobal <- system.time( ## minimally faster than the List format, but not as fast as the matrix
##     for(ii in 1:its) {
##         for(jj in 1:length(state.var.nms)) assign(paste0(state.var.nms[jj],'L'), get(state.var.nms[jj]))
##         e.sex <- sexually.activeMat[,ii]
##         s..[e.sex]   <- s..L[e.sex]*(1-p.m.bef)*(1-p.f.bef)
##         mb.a1[e.sex] <- s..L[e.sex]*p.m.bef*(1-p.f.bef)
##         mb.a2[e.sex] <- mb.a1L[e.sex]*(1 - p.f.bef)
##         mb.[e.sex]   <- mb.a2L[e.sex]*(1 - p.f.bef) + mb.L[e.sex]*(1 - p.f.bef)
##         f.ba1[e.sex] <- s..L[e.sex]*p.f.bef*(1-p.m.bef)
##         f.ba2[e.sex] <- f.ba1L[e.sex]*(1 - p.m.bef)
##         f.b[e.sex]   <- f.ba2L[e.sex]*(1 - p.m.bef) + f.bL[e.sex]*(1 - p.m.bef)
##         hb1b2[e.sex] <- hb1b2L[e.sex] + .5 * s..L[e.sex]*p.m.bef*p.f.bef + (mb.a1L[e.sex] + mb.a2L[e.sex] + mb.L[e.sex]) * p.f.bef
##         hb2b1[e.sex] <- hb2b1L[e.sex] + .5 * s..L[e.sex]*p.m.bef*p.f.bef + (f.ba1L[e.sex] + f.ba2L[e.sex] + f.bL[e.sex]) * p.m.bef
##     }
##     )
## serostatesList2 <- list(s.., mb.a1, mb.a2, mb., f.ba1, f.ba2, f.b, hb1b2, hb2b1)
## names(serostatesList2) <- state.var.nms
orig <- rbind(tDT, tMat, tMatInd, tList)#, tglobal)
#rbind(tDT/tDT, tMat/tDT, tMatInd/tDT,tList/tDT)#, tglobal/tDT)
## signif(cbind(serostates[1:10,s..],serostatesMat[1:10,'s..'], serostatesList[['s..']][1:10], serostatesList2[['s..']][1:10]),3)

pre.coupleDTc <- cmpfun(pre.coupleDT)
tDT <- system.time(
    for(ii in 1:its) {
        serostates <- pre.coupleDTc(serostates, sexually.activeMat[,ii])
    }
    ) ## about 2.25 seconds
pre.coupleMatc <- cmpfun(pre.coupleMat)
tMat <- system.time(
    for(ii in 1:its) {
        serostatesMat <- pre.coupleMatc(serostatesMat, sexually.activeMat[,ii])
    }
    ) ## about 6 seconds
pre.coupleMatIndc <- cmpfun(pre.coupleMatInd)
tMat <- system.time(
    for(ii in 1:its) {
        serostatesMatInd <- pre.coupleMatIndc(serostatesMatInd, sexually.activeMat[,ii])
    }
    ) ## about 6 seconds
pre.coupleMatOldc <- cmpfun(pre.coupleMatOld)
tList <- system.time(
    for(ii in 1:its) {
        serostatesList <- pre.coupleMatOldc(serostatesList, sexually.activeMat[,ii])
    }
    ) ## about 9 seconds
origcmp <- rbind(tDT, tMat, tMatInd, tList)#, tglobal)
#rbind(tDT/tDT, tMat/tDT, tMatInd/tDT,tList/tDT)#, tglobal/tDT)
## signif(cbind(serostates[1:10,s..],serostatesMat[1:10,'s..'], serostatesList[['s..']][1:10], serostatesList2[['s..']][1:10]),3)

orig
origcmp

Rprof(pre.coupleMat(serostatesMat, sexually.activeMat[,ii]))
