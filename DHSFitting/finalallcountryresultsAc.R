## Create Summary Figures of all the countries results
rm(list=ls())
library(abind)
library(coda)
result.dir <- "~/Dropbox/Inactive Projects/Finished Projects/Lancet MS/Couples Model Revision 121013/Output/sims/runs/"

#result.dir <- "~/Documents/R files/discordant couples/sims/runs/"
setwd(result.dir)
dirnms <- list.files()
dirnms <- dirnms[!grepl("sim",dirnms)]

#dirnms <- dirnms[!grepl("parv",dirnms) & !grepl("lcarv",dirnms) &!grepl("fsd",dirnms)] #lcarv
#dirnms <- dirnms[grepl("parv",dirnms) & !grepl("lcarv",dirnms)] #parv
#dirnms <- dirnms[grepl("parv",dirnms) & grepl("lcarv",dirnms)] #parv-lcarv
#dirnms <- dirnms[grepl("fsd",dirnms)] # fsd
result.out.dir <- "~/Dropbox/Inactive Projects/Finished Projects/Lancet MS/Couples Model Revision 121013/Output/main/"
ndir <- length(dirnms)
sigmas <- array(NA, c(6,6,10))
parnames <- c("bmb","bfb","bme","bfe","bmp","bfp","lrho")
haznames <- c("bmb","bfb","bme","bfe","bmp","bfp")
do.extra <- F

## for(ii in 1:ndir) file.rename(dirnms[ii], ~/Documents/R files/discordant couples/sims/runs old/,

for(c.ind in 1:ndir) {
    ##  do this in the loop bc loading workspaces which is dangerous
    tempdir <- paste(result.dir, dirnms[c.ind], "/files/", sep ="")
    setwd(tempdir)
    load("workspace.Rdata")
    sigmas[,,c.ind] <- sigma
    ## get beta for annual rates
    hazpars <- pars[parnames,]
    ## get geometric mean parameters for each transmission route
    bmb.vec <- unlist(mcmc.out[,"bmb"])
    bfb.vec <- unlist(mcmc.out[,"bfb"])
    bmp.vec <- unlist(mcmc.out[,"bmp"])
    bfp.vec <- unlist(mcmc.out[,"bfp"])    
    bme.vec <- unlist(mcmc.out[,"bme"])
    bfe.vec <- unlist(mcmc.out[,"bfe"])
    bb <- sqrt(bmb.vec*bfb.vec)
    be <- sqrt(bme.vec*bfe.vec)
    bp <- sqrt(bmp.vec*bfp.vec)
    rr.ep <- be / bp
    rr.bp <- bb / bp
    rr.bp.m <- bmb.vec / bp
    rr.bp.f <- bfb.vec / bp
    rr.eb <- be / bb
    gmpars <- data.frame(bb=bb,be=be,bp=bp, rr.ep = rr.ep, rr.bp = rr.bp, rr.eb = rr.eb, rr.bp.m, rr.bp.f)
    gmpars <- t(apply(gmpars, 2, function(x) quantile(x, c(.025,.5,.975))))
    pars <- rbind(gmpars, pars)      
    if(c.ind==1)
      {
        p.val.vec <- p.val
        names(p.val.vec)[length(p.val.vec)] <- as.character(dat$group[1])
        ## create paramater array for all countries
        pars.arr <- as.array(pars)
        temp.fr <- data.frame(bmb = bmb.vec, bfb = bfb.vec, bme = bme.vec, bfe = bfe.vec, bmp = bmp.vec, bfp = bfp.vec,
                              bb=bb,be=be,bp=bp, rr.ep = rr.ep, rr.bp = rr.bp, rr.eb = rr.eb, rr.bp.m, rr.bp.f)
        chains.arr <- array(temp.fr)
      }else{
        p.val.vec <- c(p.val.vec, p.val)
        names(p.val.vec)[length(p.val.vec)] <- as.character(dat$group[1])
        pars.arr <- abind(pars.arr,pars,along=3)
        temp.fr <- data.frame(bmb = bmb.vec, bfb = bfb.vec, bme = bme.vec, bfe = bfe.vec, bmp = bmp.vec, bfp = bfp.vec,
                              bb=bb,be=be,bp=bp, rr.ep = rr.ep, rr.bp = rr.bp, rr.eb = rr.eb, rr.bp.m, rr.bp.f)
        chains.arr <- abind(chains.arr, temp.fr, along = 3)
      }
    ## creade median (CI) table
    char.pars <- signif(pars,3)
    char.pars <- paste(char.pars[,2], " (", char.pars[,1], ", ", char.pars[,3],")", sep = "")
    if(c.ind==1) {
      allpars.tab <- data.frame(char.pars)
      rownames(allpars.tab) <- rownames(pars)
    }else{
      allpars.tab <- cbind(allpars.tab,char.pars)
    }
    colnames(allpars.tab)[c.ind] <- as.character(dat$group[1])
  }
setwd(result.out.dir)
save(pars.arr, file = "pars.arr.Rdata")
save(chains.arr, file = "chains.arr.Rdata")
ds.nm <- colnames(allpars.tab)
save(ds.nm, file = "ds.name.Rdata")     # save data set names (i.e. 3rd dimension index for pars.arr)
## reverse order so its in alphabetical order on barplots
pars.arr <- pars.arr[,,dim(pars.arr)[3]:1]
ds.nm <- rev(ds.nm)
ds.nm[ds.nm=="WA"] <- "West Africa"


allpars.tab <- data.frame(par = rownames(pars), allpars.tab)
save(sigmas, file="sigmas.Rdata")
write.csv(allpars.tab, file="all country results.csv")

cis <- function(x) {
  x <- signif(x,2)
  apply(x, 1, function(x) paste(x[2], "\n (", x[1], ", ", x[3], ")", sep = ""))
}


hazs <- c("bmb","bfb","bme","bfe","bmp","bfp")
pars.arr[hazs,,1]
haztab <- t(apply(12*pars.arr[hazs,,], 3, cis)) # to make them annual rates
rownames(haztab) <- ds.nm
setwd(result.out.dir)
write.csv(haztab[nrow(haztab):1,], file ="pars table.csv")

## geometric means
ghazs <- c('bb','be','bp')
pars.arr[ghazs,,1]
ghaztab <- t(apply(12*pars.arr[ghazs,,], 3, cis)) # to make them annual rates
rownames(ghaztab) <- ds.nm
setwd(result.out.dir)
write.csv(ghaztab[nrow(ghaztab):1,], file ="gpars table.csv")

rrs <- c("rr.m.out","rr.f.out","rr.m.in","rr.f.in","rr.m.pbef","rr.f.pbef", "rr.mf.exc","rr.mf.bef","rho","rr.mf.exc.cont","rr.mf.bef.cont")
pars.arr <- abind(pars.arr, rho = exp(pars.arr["lrho",,]), along = 1)
rrtab <- t(apply(pars.arr[rrs,,], 3, cis)) # to make them annual rates
rownames(rrtab) <- ds.nm
write.csv(rrtab[nrow(rrtab):1,], file ="rrs table.csv")

save(p.val.vec, file = "p.val.vec.Rdata")



###################################################################### 
pdf("serostatus breakdown.pdf", width = 5.5, height = 4)
par(mar=c(5,6,.5,.5))
load("~/Dropbox/disc mod backup/Lancet MS/Couples Model Revision 121013/R scripts & Data Files/alldhs.Rdata")
tab1 <- xtabs(~group + ser, dat)
tab1 <- tab1[nrow(tab1):1,]
total <- rowSums(tab1)
total <- as.matrix(total)[,rep(1,4)]
tab1 <- tab1/total
tab1
rownames(tab1)[rownames(tab1)=="WA"] <- "West Africa"
cols <- c("yellow","green","purple","dark gray")
barplot(t(tab1), beside = F, names.arg = rownames(tab1), horiz = T, las = 2,
        col = cols, xlab = "proportion of couples in serogroup", border = NA)
dev.off()
###################################################################### 
tab1 <- tab1[nrow(tab1):1,]
write.csv(signif(tab1,3)*100, "serostatus breakdown.csv")

pdf("serostatus breakdown leg.pdf", width = 6, height = 2.5)
par(mar=rep(0,4))
plot(0,0, type = "n", bty = "n", xlab = "", ylab = "", axes = F)
legend("top", c("M+ F+", "M+ F-", "M- F+", "M- F-"),
       pch = 15, bty = "n", col = cols, cex = 1.5)
dev.off()

###################################################################### 
pdf("serostatus breakdown neg only.pdf", width = 8, height = 5)
par(mar=c(4,6,3,2))
tab1 <- xtabs(~group + ser, dat)
tab1 <- tab1[nrow(tab1):1,]
total <- rowSums(tab1)
total <- as.matrix(total)[,rep(1,4)]
tab1 <- tab1/total
tab1
rownames(tab1)[rownames(tab1)=="WA"] <- "West Africa"
cols <- c(NA,"green","purple","dark gray")
barplot(t(tab1), beside = F, names.arg = rownames(tab1), horiz = T, las = 2,
        col = cols, xlab = "proportion of couples in serogroup", border = NA)
dev.off()
###################################################################### 

 
###################################################################### 
## beta's by country
for(ff in 1:3)
  {
    pdf(paste('beta before',ff,'.pdf'), width = 8, height = 6)
    ff <- 1:ff
    xlim <- c(0, 12*max(pars.arr[hazs,,]))
    par(mar=c(5,8,2,2))
    plot(0,0, type = "n", xlim = xlim, xlab = expression(beta[b]),
         bty = "n", ylim = c(.4, dim(pars.arr)[3]+.2), yaxt="n", ylab = "",
         main="")
    if(1 %in% ff)
      {
        points(12*pars.arr["bmb",2,], 1:dim(pars.arr)[3], pch = 15)
        arrows(12*pars.arr["bmb",1,], 1:dim(pars.arr)[3],
               12*pars.arr["bmb",3,], 1:dim(pars.arr)[3], 
               angle = 90, length=.05, code = 3, lwd = 2)
      }
    if(2 %in% ff)
      {
        points(12*pars.arr["bfb",2,], 1:dim(pars.arr)[3]-.3 , pch = 15,
               col = "brown")
        arrows(12*pars.arr["bfb",1,], 1:dim(pars.arr)[3]-.3,
               12*pars.arr["bfb",3,], 1:dim(pars.arr)[3]-.3, 
               angle = 90, length=.05, code = 3, lwd = 2, col = "brown")
      }
    if(3 %in% ff)
      {
        points(12*pars.arr["bb",2,], 1:dim(pars.arr)[3] - .5, pch = 15, col = 'orange')
        arrows(12*pars.arr["bb",1,], 1:dim(pars.arr)[3] - .5, 
               12*pars.arr["bb",3,], 1:dim(pars.arr)[3] - .5, col = 'orange',
               angle = 90, length=.05, code = 3, lwd = 2)
      }
    axis(2, at = 1:dim(pars.arr)[3]+.15, lab = ds.nm, las = 2)
    if(show.geom)
      {
        legend("topright", c("male","female",'geometric mean'), pch = 15, col = c("black","brown",'orange'), bty = "n")
      }else{
        legend("topright", c("male","female"), pch = 15, col = c("black","brown"), bty = "n")
      }
    dev.off()
    pdf("beta extracouple.pdf", width = 8, height = 6)
    par(mar=c(5,8,2,2))
    plot(0,0, type = "n", xlim = xlim, xlab = expression(beta[e]),
         bty = "n", ylim = c(.4, dim(pars.arr)[3]+.2), yaxt="n", ylab = "",
         main="")
    if(1 %in% ff)
      {
        points(12*pars.arr["bme",2,], 1:dim(pars.arr)[3], pch = 15)
        arrows(12*pars.arr["bme",1,], 1:dim(pars.arr)[3],
               12*pars.arr["bme",3,], 1:dim(pars.arr)[3], 
               angle = 90, length=.05, code = 3, lwd = 2)
      }
    if(2 %in% ff)
      {
        points(12*pars.arr["bfe",2,], 1:dim(pars.arr)[3]-.3 , pch = 15,
               col = "brown")
        arrows(12*pars.arr["bfe",1,], 1:dim(pars.arr)[3]-.3,
               12*pars.arr["bfe",3,], 1:dim(pars.arr)[3]-.3, 
               angle = 90, length=.05, code = 3, lwd = 2, col = "brown")
      }
    if(3 %in% ff)
      {
        points(12*pars.arr["be",2,], 1:dim(pars.arr)[3] - .5, pch = 15, col = 'orange')
        arrows(12*pars.arr["be",1,], 1:dim(pars.arr)[3] - .5, 
               12*pars.arr["be",3,], 1:dim(pars.arr)[3] - .5, col = 'orange',
               angle = 90, length=.05, code = 3, lwd = 2)
        legend("topright", c("male","female",'geometric mean'), pch = 15, col = c("black","brown",'orange'), bty = "n")
      }else{
                legend("topright", c("male","female"), pch = 15, col = c("black","brown"), bty = "n")
              }
    axis(2, at = 1:dim(pars.arr)[3]+.15, lab = ds.nm, las = 2)
    dev.off()
    pdf("beta partner.pdf", width = 8, height = 6)
    par(mar=c(5,8,2,2))
    plot(0,0, type = "n", xlim = xlim, xlab = expression(beta[p]),
         bty = "n", ylim = c(.4, dim(pars.arr)[3]+.2), yaxt="n", ylab = "",
         main="")
    points(12*pars.arr["bmp",2,], 1:dim(pars.arr)[3], pch = 15)
    arrows(12*pars.arr["bmp",1,], 1:dim(pars.arr)[3],
           12*pars.arr["bmp",3,], 1:dim(pars.arr)[3], 
           angle = 90, length=.05, code = 3, lwd = 2)
    points(12*pars.arr["bfp",2,], 1:dim(pars.arr)[3]-.3 , pch = 15,
           col = "brown")
    arrows(12*pars.arr["bfp",1,], 1:dim(pars.arr)[3]-.3,
           12*pars.arr["bfp",3,], 1:dim(pars.arr)[3]-.3, 
           angle = 90, length=.05, code = 3, lwd = 2, col = "brown")
    axis(2, at = 1:dim(pars.arr)[3]+.15, lab = ds.nm, las = 2)
    if(show.geom)
      {
        points(12*pars.arr["bp",2,], 1:dim(pars.arr)[3] - .5, pch = 15, col = 'orange')
        arrows(12*pars.arr["bp",1,], 1:dim(pars.arr)[3] - .5, 
               12*pars.arr["bp",3,], 1:dim(pars.arr)[3] - .5, col = 'orange',
               angle = 90, length=.05, code = 3, lwd = 2)
      }
    if(show.geom)
      {
        legend("topright", c("male","female",'geometric mean'), pch = 15, col = c("black","brown",'orange'), bty = "n")
      }else{
        legend("topright", c("male","female"), pch = 15, col = c("black","brown"), bty = "n")
      }
    dev.off()
###################################################################### 
  }

###################################################################### 
## geometric mean beta's by country
xlim <- c(0, 12*max(pars.arr[hazs,,]))
pdf("gbeta before.pdf", width = 8, height = 6)
par(mar=c(5,8,2,2))
plot(0,0, type = "n", xlim = xlim, xlab = expression(beta[b]),
     bty = "n", ylim = c(.4, dim(pars.arr)[3]+.2), yaxt="n", ylab = "",
     main="")
points(12*pars.arr["bb",2,], 1:dim(pars.arr)[3], pch = 15)
arrows(12*pars.arr["bb",1,], 1:dim(pars.arr)[3],
       12*pars.arr["bb",3,], 1:dim(pars.arr)[3], 
       angle = 90, length=.05, code = 3, lwd = 2)
axis(2, at = 1:dim(pars.arr)[3]+.15, lab = ds.nm, las = 2)
dev.off()
#xlim <- c(0,.4)
pdf("gbeta extracouple.pdf", width = 6, height = 6)
par(mar=c(5,8,2,2))
plot(0,0, type = "n", xlim = xlim, xlab = expression(beta[e]),
     bty = "n", ylim = c(.4, dim(pars.arr)[3]+.2), yaxt="n", ylab = "",
     main="")
points(12*pars.arr["be",2,], 1:dim(pars.arr)[3], pch = 15)
arrows(12*pars.arr["be",1,], 1:dim(pars.arr)[3],
       12*pars.arr["be",3,], 1:dim(pars.arr)[3], 
       angle = 90, length=.05, code = 3, lwd = 2)
axis(2, at = 1:dim(pars.arr)[3]+.15, lab = ds.nm, las = 2)
dev.off()
pdf("gbeta partner.pdf", width = 6, height = 6)
par(mar=c(5,8,2,2))
plot(0,0, type = "n", xlim = xlim, xlab = expression(beta[p]),
     bty = "n", ylim = c(.4, dim(pars.arr)[3]+.2), yaxt="n", ylab = "",
     main="")
points(12*pars.arr["bp",2,], 1:dim(pars.arr)[3], pch = 15)
arrows(12*pars.arr["bp",1,], 1:dim(pars.arr)[3],
       12*pars.arr["bp",3,], 1:dim(pars.arr)[3], 
       angle = 90, length=.05, code = 3, lwd = 2)
axis(2, at = 1:dim(pars.arr)[3]+.15, lab = ds.nm, las = 2)
dev.off()
###################################################################### 


###################################################################### 
## beta ratios by country
## 
pdf("beta ratio extracouple to before.pdf", width = 6, height = 6)
par(mar=c(5,8,2,2))
plot(1,0, type = "n", xlim = c(10^-2,10^2), xlab = expression(12*beta[before]),
     bty = "n", ylim = c(.4, dim(pars.arr)[3]+.2), yaxt="n", ylab = "",
     log="x",
     main="transmissibility ratio: extracouple / before couple formation")
points(pars.arr["rr.m.out",2,], 1:dim(pars.arr)[3], pch = 15)
arrows(pars.arr["rr.m.out",1,], 1:dim(pars.arr)[3],
       pars.arr["rr.m.out",3,], 1:dim(pars.arr)[3], 
       angle = 90, length=.05, code = 3, lwd = 2)
points(pars.arr["rr.f.out",2,], 1:dim(pars.arr)[3]-.3 , pch = 15,
       col = "brown")
arrows(pars.arr["rr.f.out",1,], 1:dim(pars.arr)[3]-.3,
       pars.arr["rr.f.out",3,], 1:dim(pars.arr)[3]-.3, 
       angle = 90, length=.05, code = 3, lwd = 2, col = "brown")
axis(2, at = 1:dim(pars.arr)[3]+.15, lab = ds.nm, las = 2)
legend("topleft", c("male","female"), pch = 15, col = c("black","brown"),
       bty = "n")
axis(1, at = c(10^c(-4:-2),1/1:10,2:10,10^c(2:4)), lab = NA)
abline(v=1, lty = 2)
dev.off()


###################################################################### 
pdf("beta ratio m to f extracouple.pdf", width = 6, height = 6)
par(mar=c(5,8,2,2))
plot(1,0, type = "n", xlim = c(10^-2,10^2), xlab = expression(12*beta[before]),
     bty = "n", ylim = c(.4, dim(pars.arr)[3]+.2), yaxt="n", ylab = "",
     log="x",
     main="extracouple transmissibility ratio: m:f")
points(pars.arr["rr.mf.exc",2,], 1:dim(pars.arr)[3], pch = 15)
arrows(pars.arr["rr.mf.exc",1,], 1:dim(pars.arr)[3],
       pars.arr["rr.mf.exc",3,], 1:dim(pars.arr)[3], 
       angle = 90, length=.05, code = 3, lwd = 2)
axis(2, at = 1:dim(pars.arr)[3]+.15, lab = ds.nm, las = 2)
axis(1, at = c(10^c(-4:-2),1/1:10,2:10,10^c(2:4)), lab = NA)
abline(v=1, lty = 2)
dev.off()
###################################################################### 
pdf("beta ratio m to f extracouple xRHO.pdf", width = 6, height = 6)
par(mar=c(5,8,2,2))
plot(1,0, type = "n", xlim = c(10^-2,10^2), xlab = expression(12*beta[before]),
     bty = "n", ylim = c(.4, dim(pars.arr)[3]+.2), yaxt="n", ylab = "",
     log="x",
     main="extracouple transmissibility ratio: m:f")
points(pars.arr["rr.mf.exc.cont",2,], 1:dim(pars.arr)[3], pch = 15)
arrows(pars.arr["rr.mf.exc.cont",1,], 1:dim(pars.arr)[3],
       pars.arr["rr.mf.exc.cont",3,], 1:dim(pars.arr)[3], 
       angle = 90, length=.05, code = 3, lwd = 2)
axis(2, at = 1:dim(pars.arr)[3]+.15, lab = ds.nm, las = 2)
axis(1, at = c(10^c(-4:-2),1/1:10,2:10,10^c(2:4)), lab = NA)
abline(v=1, lty = 2)
dev.off()
######################################################################

###################################################################### 
pdf("beta ratio m to f before.pdf", width = 6, height = 6)
par(mar=c(5,8,2,2))
plot(1,0, type = "n", xlim = c(10^-2,10^2), xlab = expression(12*beta[before]),
     bty = "n", ylim = c(.4, dim(pars.arr)[3]+.2), yaxt="n", ylab = "",
     log="x",
     main="before couple formation transmissibility ratio: m:f")
points(pars.arr["rr.mf.bef",2,], 1:dim(pars.arr)[3], pch = 15)
arrows(pars.arr["rr.mf.bef",1,], 1:dim(pars.arr)[3],
       pars.arr["rr.mf.bef",3,], 1:dim(pars.arr)[3], 
       angle = 90, length=.05, code = 3, lwd = 2)
axis(2, at = 1:dim(pars.arr)[3]+.15, lab = ds.nm, las = 2)
axis(1, at = c(10^c(-4:-2),1/1:10,2:10,10^c(2:4)), lab = NA)
abline(v=1, lty = 2)
dev.off()
###################################################################### 
pdf("beta ratio m to f before xRHO.pdf", width = 6, height = 6)
par(mar=c(5,8,2,2))
plot(1,0, type = "n", xlim = c(10^-2,10^2), xlab = expression(12*beta[before]),
     bty = "n", ylim = c(.4, dim(pars.arr)[3]+.2), yaxt="n", ylab = "",
     log="x",
     main="before couple formation transmissibility ratio: m:f")
points(pars.arr["rr.mf.bef.cont",2,], 1:dim(pars.arr)[3], pch = 15)
arrows(pars.arr["rr.mf.bef.cont",1,], 1:dim(pars.arr)[3],
       pars.arr["rr.mf.bef.cont",3,], 1:dim(pars.arr)[3], 
       angle = 90, length=.05, code = 3, lwd = 2)
axis(2, at = 1:dim(pars.arr)[3]+.15, lab = ds.nm, las = 2)
axis(1, at = c(10^c(-4:-2),1/1:10,2:10,10^c(2:4)), lab = NA)
abline(v=1, lty = 2)
dev.off()
######################################################################
###################################################################### 


###################################################################### 
pdf("Figure 4 - barplot prob extracouple 12 months.pdf", width = 4, height = 3)
cex.names <- .7
par(mar=c(3,.8,1,.5), mfrow = c(1,2), oma = c(0,3.5,0,0))
cols <- c("red") #,"dodgerblue1")
bp.temp <- barplot(pars.arr["prop.exc.m",2,], beside = F, xlim = c(0,1), names.arg = ds.nm, horiz=2, border = NA,
                   las = 2, bty = "n", ylab = "", main="", col = cols, cex.names = cex.names, cex.axis = cex.names)
mtext("males", line = 0, side = 3, cex = cex.names) 
arrows(pars.arr["prop.exc.m",1,], bp.temp,
       pars.arr["prop.exc.m",3,], bp.temp,
       angle = 90, length=.05, code = 3, lwd = 1)
bp.temp <- barplot(pars.arr["prop.exc.f",2,], xlim = c(0,1), names.arg = rep("",length(ds.nm)), horiz=2, border = NA,
                   las = 2, bty = "n", ylab = "", main="", col = cols, cex.axis = cex.names)
mtext("females", line = 0, side = 3, cex = cex.names) 
arrows(pars.arr["prop.exc.f",1,], bp.temp,
       pars.arr["prop.exc.f",3,], bp.temp,
       angle = 90, length=.05, code = 3, lwd = 1)
mtext("proportion of HIV transmission in next year due to extra-couple sex", side = 1, line = -1.2,
      outer = T, at = .4, cex = cex.names)
dev.off()


###################################################################### 
######################################################################
## future transmission break-down multiple bars
###################################################################### 
## to have barplot be in alphabetical order
load("~/Dropbox/disc mod backup/Lancet MS/Couples Model Revision 121013/R scripts & Data Files/alldhs.Rdata")
ds.nm.wa <- ds.nm
## model calculates incidence per 1000 couples where either partner is at risk. Will change pcalc.R later but for now adjust to incidence per 1000 negative individuals
tsamptab <- xtabs(~group, dat, subset = ser!=1)
names(tsamptab)[names(tsamptab)=="WA"] <- "West Africa"
tnmatrix <- array(rep(tsamptab[match(ds.nm.wa, names(tsamptab))], each = 12),
                  dim = c(4,3, length(ds.nm.wa)))
## # of couples with negative mles
show <- c("n.m.exc.cc", "n.m.exc.dc", "n.m.part.dc","n.m.part.cc")
msamptab <- xtabs(~group, dat, subset = ser %in% c(4,3))
names(msamptab)[names(msamptab)=="WA"] <- "West Africa"
mnmatrix <- array(rep(msamptab[match(ds.nm.wa, names(msamptab))], each = 12),
                  dim = c(4,3, length(ds.nm.wa)))
mbdn <- pars.arr[show,,]*tnmatrix/mnmatrix
show <- c("n.m.exc.tot", "n.m.part.tot")
mbdn.tot <- pars.arr[show,,]*tnmatrix[1:2,,]/mnmatrix[1:2,,]
## # of couples with negative mles
show <- c("n.f.exc.cc", "n.f.exc.dc", "n.f.part.dc","n.f.part.cc")
fsamptab <- xtabs(~group, dat, subset = ser %in% c(4,3))
names(fsamptab)[names(fsamptab)=="WA"] <- "West Africa"
fnmatrix <- array(rep(fsamptab[match(ds.nm.wa, names(fsamptab))], each = 12),
                  dim = c(4,3, length(ds.nm.wa)))
fbdn <- pars.arr[show,,]*tnmatrix/fnmatrix
show <- c("n.f.exc.tot", "n.f.part.tot")
fbdn.tot <- pars.arr[show,,]*tnmatrix[1:2,,]/fnmatrix[1:2,,]
## give incidence/100 by dividing by sample size
mbdn.fut <- mbdn/10                         # per 100 person years instead of 1000
fbdn.fut <- fbdn/10                         # per 100 person years instead of 1000
mbdn.fut.tot <- mbdn.tot/10                         # per 100 person years instead of 1000
fbdn.fut.tot <- fbdn.tot/10                         # per 100 person years instead of 1000


######################################################################
## dc+ cc breakdown in Observed
beside <- F
show <- c("pibUA","pieUA","pipUA")
mbdn.all <- pars.arr[show,,]
mbdn.all[,2,] <- apply(mbdn.all[,2,], 2, function(x) {x/sum(x)}) # normalize medians to add to 1 (don't have to b/c they come from accross chains)
show <- c("piUbA","piUeA","piUpA")
fbdn.all <- pars.arr[show,,]
fbdn.all[,2,] <- apply(fbdn.all[,2,], 2, function(x) {x/sum(x)})
space <- .2
if(beside) space <- c(0,.5)
###################################################################### 
## plot breakdown of all OBSERVED Infections
pdf("cc & dc breakdown.pdf", width = 8, height = 5)
par(mar=c(4,6,3,2), mfrow = c(1,2), oma = c(0,0,1,0))
cols <- c("dark gray","red", "blue")
xlim <- c(0,1)
bp.temp <- barplot(mbdn.all[,2,], beside = beside, names.arg = ds.nm.wa, las =2, space = space,
                   horiz = T, col = cols, xlim = xlim, main = "males")
if(beside) arrows(mbdn.all[,1,], bp.temp, mbdn.all[,3,], bp.temp,
                  angle = 90, length=.025, code = 3, lwd = 2)
bp.temp <- barplot(fbdn.all[,2,], beside = beside, names.arg = ds.nm.wa, las =2, space = space,
                   horiz = T, col = cols, xlim = xlim, main = "females")
if(beside) arrows(fbdn.all[,1,], bp.temp, fbdn.all[,3,], bp.temp,
                  angle = 90, length=.025, code = 3, lwd = 2)
mtext("proportion of infected individuals", side = 1, line = -1.3, outer=T,
      cex = 1.3)
## if(survive)
##   {
##     mtext("with AIDS deaths", side = 3, line = -.5, outer = T, cex = 1.3)
##   }else{
##     mtext("without AIDS deaths", side = 3, line = -.5, outer = T, cex = 1.3)
##   }
                                        #legend("bottomright", c("before", "extracouply", "infected by partner"), pch = 15, bty = "n", col = cols, cex = .7)
dev.off()



######################################################################
## concordant + breakdown in OBSERVED couples
beside <- F
pdf("CC breakdown.pdf", width = 8, height = 5)
par(mar=c(4,6,3,2), mfrow = c(1,2), oma = c(0,0,1,0))
cols <- c("dark gray","red","blue")
xlim <- c(0,1)
show <- c("pibPA", "piePA", "pipPA")
mbdn.cc <- pars.arr[show,,]
mbdn.cc[,2,] <- apply(mbdn.cc[,2,], 2, function(x) {x/sum(x)})
show <- c("piPbA", "piPeA", "piPpA")
fbdn.cc <- pars.arr[show,,]
fbdn.cc[,2,] <- apply(fbdn.cc[,2,], 2, function(x) {x/sum(x)})
space <- .2
if(beside) space <- c(0,.5)
bp.temp <- barplot(mbdn.cc[,2,], beside = beside, names.arg = ds.nm.wa, las =2, space = space,
                   horiz = T, col = cols, xlim = xlim, main = "males")
if(beside) arrows(mbdn.cc[,1,], bp.temp, mbdn.cc[,3,], bp.temp,
                  angle = 90, length=.025, code = 3, lwd = 2)
bp.temp <- barplot(fbdn.cc[,2,], beside = beside, names.arg = ds.nm.wa, las =2, space = space,
                   horiz = T, col = cols, xlim = xlim, main = "females")
if(beside) arrows(fbdn.cc[,1,], bp.temp, fbdn.cc[,3,], bp.temp,
                  angle = 90, length=.025, code = 3, lwd = 2)
mtext("proportion of concordant positive couples", side = 1, line = -1.3, outer=T,
      cex = 1.3)
## if(survive)
##   {
##     mtext("with AIDS deaths", side = 3, line = -.5, outer = T, cex = 1.3)
##   }else{
##     mtext("without AIDS deaths", side = 3, line = -.5, outer = T, cex = 1.3)
##   }
dev.off()
pdf("cc breakdown leg.pdf", width = 6, height = 2.5)
par(mar=rep(0,4))
plot(0,0, type = "n", bty = "n", xlab = "", ylab = "", axes = F)
legend("top", c("infected before couple formation", "infected extracouply",
                "infected by partner"),
       pch = 15, bty = "n", col = cols, cex = 1.5)
dev.off()
######################################################################

######################################################################
## disconcordant + breakdown in OBSERVED COUPLES
beside <- F
pdf("dc breakdown.pdf", width = 8, height = 5)
par(mar=c(4,6,3,2), mfrow = c(1,2), oma = c(0,0,1,0))
show <- c("pibNA","pieNA")
mbdn.dc <- pars.arr[show,,]
mbdn.dc[,2,] <- apply(mbdn.dc[,2,], 2, function(x) {x/sum(x)})
show <- c("piNbA","piNeA")
fbdn.dc <- pars.arr[show,,]
fbdn.dc[,2,] <- apply(fbdn.dc[,2,], 2, function(x) {x/sum(x)})
space <- .2
if(beside) space <- c(0,.5)
## plot
cols <- c("dark gray","red")
xlim <- c(0,1)
bp.temp <- barplot(mbdn.dc[,2,], beside = beside, names.arg = ds.nm.wa, las =2, space = space,
                   horiz = T, col = cols, xlim = xlim, main = "males")
if(beside) arrows(mbdn.dc[,1,], bp.temp, mbdn.dc[,3,], bp.temp,
                  angle = 90, length=.025, code = 3, lwd = 2)
bp.temp <- barplot(fbdn.dc[,2,], beside = beside, names.arg = ds.nm.wa, las =2, space = space,
                   horiz = T, col = cols, xlim = xlim, main = "females")
if(beside) arrows(fbdn.dc[,1,], bp.temp, fbdn.dc[,3,], bp.temp,
                  angle = 90, length=.025, code = 3, lwd = 2)
mtext("proportion of discordant positive couples", side = 1, line = -1.3, outer=T,
      cex = 1.3)
## if(survive)
##   {
##     mtext("with AIDS deaths", side = 3, line = -.5, outer = T, cex = 1.3)
##   }else{
##     mtext("without AIDS deaths", side = 3, line = -.5, outer = T, cex = 1.3)
##   }
dev.off()
pdf("dc breakdown leg.pdf", width = 6, height = 2.5)
par(mar=rep(0,4))
plot(0,0, type = "n", bty = "n", xlab = "", ylab = "", axes = F)
legend("top", c("infected before couple formation", "infected extracouply"),
       pch = 15, bty = "n", col = cols, cex = 1.5)
dev.off()
###################################################################### 

######################################################################
## Index partner breakdown, gender pooled, estimator sumI
######################################################################
pdf("disc breakdown gender pooled estimator1.pdf", width = 4.21, height = 2.5)
par(mar=c(4.5,3.5,1.5,.5))
show <- c("piGb1.sumI","piGe1.sumI","piG.b1sumI","piG.e1sumI")
bdn.index <- pars.arr[show,,]
bdn.index[,2,] <- apply(bdn.index[,2,], 2, function(x) {x/sum(x)})
xlim <- c(0,1)
space <- 0
beside <- F
if(beside) space <- c(0,.5)
cols <- c("dark gray","red", NA, NA)
bp.temp <- barplot(bdn.index[,2,], beside = beside, names.arg = ds.nm.wa, las =2, space = space,
                   horiz = T, col = cols, xlim = xlim, main = "males")
par(new=T)
cols <- c("dark gray", "red","dark gray","red")
bp.temp <- barplot(bdn.index[,2,], beside = beside, names.arg = rep(NA, length(ds.nm.wa)), las =2, space = space,
                   horiz = T, col = cols, xlim = xlim, main = "males", density = 30, axes = F)
if(beside) arrows(bdn.index[,1,], bp.temp, bdn.index[,3,], bp.temp,
                  angle = 90, length=.025, code = 3, lwd = 2)
mtext("proportion of index infections in couples", side = 1, line = -1.3, outer=T,
      cex = 1.3)
dev.off()

## ######################################################################
## ## Index partner breakdown, gender pooled Estimator sumIS
## ######################################################################
## pdf("disc breakdown gender pooled estimator2.pdf", width = 4.21, height = 2.5)
## par(mar=c(4.5,3.5,1.5,.5))
## show <- c("piGb1.sumIS","piGe1.sumIS","piG.b1sumIS","piG.e1sumIS")
## bdn.index <- pars.arr[show,,]
## bdn.index[,2,] <- apply(bdn.index[,2,], 2, function(x) {x/sum(x)})
## xlim <- c(0,1)
## space <- 0
## beside <- F
## if(beside) space <- c(0,.5)
## cols <- c("dark gray","red", NA, NA)
## bp.temp <- barplot(bdn.index[,2,], beside = beside, names.arg = ds.nm.wa, las =2, space = space,
##         horiz = T, col = cols, xlim = xlim, main = "males")
## par(new=T)
## cols <- c("dark gray", "red","dark gray","red")
## bp.temp <- barplot(bdn.index[,2,], beside = beside, names.arg = rep(NA, length(ds.nm.wa)), las =2, space = space,
##         horiz = T, col = cols, xlim = xlim, main = "males", density = 30, axes = F)
## if(beside) arrows(bdn.index[,1,], bp.temp, bdn.index[,3,], bp.temp,
##        angle = 90, length=.025, code = 3, lwd = 2)
## mtext("proportion of index infections in couples", side = 1, line = -1.3, outer=T,
##       cex = 1.3)
## dev.off()

######################################################################
## Index partner breakdown, gender pooled sumI
######################################################################
pdf("Figure S6 -disc breakdown gender pooled.pdf", width = 4.21, height = 3.5)
layout(matrix(c(1,2),nr=2,nc=1), w=c(1,1),h=c(.2,.9))
par(mar=rep(0,4))
border <- "white"
plot(0,0, type = "n", axes = F, bty = "n")
cols.raw <- c("dark gray","red")
legend("topleft", leg = c("male before couple", "male extra-couple"),
       fill = cols.raw, bty = "n", border = border)
legend("topright", leg = c("female before couple", "female extra-couple"),
       fill = cols.raw, bty = "n", border = border)
legend("topright", leg = c("female before couple", "female extra-couple"),
       fill = "black", bty = "n", density = 30, border = border)
par(mar=c(4.5,5.5,0,.5))
show <- c("piGb1.sumI","piGe1.sumI","piG.b1sumI","piG.e1sumI")
bdn.index <- pars.arr[show,,]
bdn.index[,2,] <- apply(bdn.index[,2,], 2, function(x) {x/sum(x)})
xlim <- c(0,1)
space <- .1
beside <- F
if(beside) space <- c(0,.5)
cols <- c(cols.raw, cols.raw)
bp.temp <- barplot(bdn.index[,2,], beside = beside, names.arg = ds.nm.wa, las =2, space = space, border = border,
                   horiz = T, col = cols, xlim = xlim, main = "", axes = F)
par(new=T)
cols <- c(cols.raw,"black","black")
bp.temp <- barplot(bdn.index[,2,], beside = beside, names.arg = rep(NA, length(ds.nm.wa)), las =2, space = space, border = border,
                   horiz = T, col = cols, xlim = xlim, main = "", density = 30, axes = F)
if(beside) arrows(bdn.index[,1,], bp.temp, bdn.index[,3,], bp.temp,
                  angle = 90, length=.025, code = 3, lwd = 2)
m.prop <- tab1[,2]/rowSums(tab1[,2:3])
axis(1, at =c(0,.25,.5,.75,1), las = 2)
## points(m.prop[names(m.prop) %in% ds.nm], bp.temp, pch = 15, cex = .7)
mtext("proportion of index infections in couples", side = 1, line = -1.3, outer=T,
      cex = 1)
dev.off()


## for text: route-contribution ranges within couple type 
apply(mbdn.dc[,2,], 1, range)
apply(fbdn.dc[,2,], 1, range)
apply(mbdn.cc[,2,], 1, range)
apply(fbdn.cc[,2,], 1, range)
apply(mbdn.all[,2,], 1, range)
apply(fbdn.all[,2,], 1, range)
range(pars.arr["prop.exc.m",2,])
range(pars.arr["prop.exc.f",2,])

## excluding drc
apply(mbdn.dc[,2,-10], 1, range)
apply(fbdn.dc[,2,-10], 1, range)
apply(mbdn.cc[,2,-10], 1, range)
apply(fbdn.cc[,2,-10], 1, range)
apply(mbdn.all[,2,-10], 1, range)
apply(fbdn.all[,2,-10], 1, range)
apply(bdn.index[,2,-10], 1, range)
range(pars.arr["prop.exc.m",2,-10])
range(pars.arr["prop.exc.f",2,-10])
## for text: ranges of hazards
apply(pars.arr[hazs,2,-10]*12,1, range)
within <- pars.arr[c("bmp","bfp"),2,]*12
colnames(within) <- ds.nm.wa

## observed transmission breakdown table
mbdn.dc3 <- abind(mbdn.dc, array(0, dim = c(1,3,dim(mbdn.dc)[3])), along = 1)
fbdn.dc3 <- abind(fbdn.dc, array(0, dim = c(1,3,dim(fbdn.dc)[3])), along = 1)
dctab <- cbind(t(apply(mbdn.dc3, 3, cis)), t(apply(fbdn.dc3, 3, cis)))
dctab[,c(3,6)] <- "-"
rownames(dctab) <- ds.nm.wa
cctab <- cbind(t(apply(mbdn.cc, 3, cis)), t(apply(fbdn.cc, 3, cis)))
rownames(cctab) <- ds.nm.wa
bothtab <- cbind(t(apply(mbdn.all, 3, cis)), t(apply(fbdn.all, 3, cis)))
rownames(bothtab) <- ds.nm.wa
alltab <- rbind(bothtab,cctab,dctab)
alltab <- alltab[nrow(alltab):1,]
write.csv(alltab, file="all historical breakdown.csv")

indtab <- t(apply(bdn.index, 3, cis))
rownames(indtab) <- ds.nm
write.csv(indtab[nrow(indtab):1,], file = "index infection table.csv")

## proj transm table
futtab <- cbind(t(apply(mbdn.fut.tot, 3, cis)), t(apply(fbdn.fut.tot, 3, cis)))
rownames(futtab) <- ds.nm.wa
futtab <- futtab[nrow(futtab):1,]
write.csv(futtab, file="incidence future breakdown.csv")


######################################################################
## Calculate proportion of future transmission in observed discordant
## couples that is exc
## mfdc <- t(pars.arr["n.m.exc.dc",,] / (pars.arr["n.m.exc.dc",,] + pars.arr["n.m.part.dc",,]))
## ffdc <- t(pars.arr["n.f.exc.dc",,] / (pars.arr["n.f.exc.dc",,] + pars.arr["n.f.part.dc",,]))
## rownames(mfdc) <- ds.nm.wa
## rownames(ffdc) <- ds.nm.wa
## mfdc <- signif(mfdc,2)
## ffdc <- signif(ffdc,2)
## mfdc.tab <- paste(mfdc[,2], "\n (", mfdc[,1], ", ", mfdc[,3],")", sep = "")
## ffdc.tab <- paste(ffdc[,2], "\n (", ffdc[,1], ", ", ffdc[,3],")", sep = "")
## fdc.tab <- cbind(mfdc.tab, ffdc.tab)
## colnames(fdc.tab) <- c("male","female")
prtab.dc <- t(apply(pars.arr[c("prop.exc.m.dc","prop.exc.f.dc"),,],3,cis))
rownames(prtab.dc) <- ds.nm.wa
prtab.dc <- prtab.dc[nrow(futtab):1,]
write.csv(prtab.dc, file="future breakdown discordant.csv")

## proj transm table
prtab <- t(apply(pars.arr[c("prop.exc.m","prop.exc.f"),,],3,cis))
rownames(prtab) <- ds.nm.wa
prtab <- prtab[nrow(futtab):1,]
write.csv(prtab, file="proportion incidence future breakdown.csv")



######################################################################
## Table 1
## data set grouping
ind <- !duplicated(dat$ds)
tab9 <- dat[ind,c("group","ds")]
for(ii in 1:length(unique(tab9$group)))
  {
    if(ii==1)
      {
        tab10 <- data.frame(group = unique(tab9$group)[ii],
                            data.sets = paste(tab9$ds[tab9$group==unique(tab9$group)[ii]], collapse = "; "))
      }else{
        tab10 <- rbind(tab10,
                       data.frame(group = unique(tab9$group)[ii],
                                  data.sets = paste(tab9$ds[tab9$group==unique(tab9$group)[ii]], collapse = "; ")))
      }
  }
write.csv(tab10, file="ds grouping.csv")


######################################################################
## Figure 3: proportional contibution of each route for OBSERVED
## Couples for all positive & predicted proportion in future.
###################################################################### 
## plot
ords <- 1:3
space <- 0.1
len <- .01
cex.names <- .8
pdf("Figure 3 - transmission breakdown.pdf", width = 4.21, height = 4.5)
par(oma = c(1,4.5,.5,.5))
par(mar=c(2.5,.5,1.5,1.5), mfrow = c(2,2))
cols <- c("dark gray","red", "dodgerblue1")[ords] # rev because reversing color order
xlim <- c(0,1)
## all
bp.temp <- barplot(mbdn.all[ords,2,], beside = beside, names.arg = ds.nm.wa, las =2, space = space, cex.names = cex.names,
                   horiz = T, col = cols, xlim = xlim,  border = NA, axes = F)
mtext("A", side = 3, line = .2, cex = cex.names)
axis(1, seq(0,1,by=.2), las = 2, cex.axis = cex.names)
bp.temp <- barplot(fbdn.all[ords,2,], beside = beside,
                   names.arg = rep("", length(ds.nm.wa)), las =2, space = space,
                   horiz = T, col = cols, xlim = xlim,  border = NA, axes = F)
axis(1, seq(0,1,by=.2), las = 2, , cex.axis = cex.names)
mtext("B", side = 3, line = .2, cex = cex.names)
if(beside) arrows(fbdn.all[ords,1,], bp.temp, fbdn.all[ords,3,], bp.temp, cex.names = cex.names,
                  angle = 90, length=len, code = 3, lwd = 1)
###################################################################### 
cols <- c("red") #,"dodgerblue1")
bp.temp <- barplot(pars.arr["prop.exc.m",2,], beside = F, xlim = c(0,1), names.arg = ds.nm, horiz=2, border = NA,
                   las = 2, bty = "n", ylab = "", main="", col = cols, cex.names = cex.names, cex.axis = cex.names)
mtext("C", side = 3, line = .2, cex = cex.names)
arrows(pars.arr["prop.exc.m",1,], bp.temp,
       pars.arr["prop.exc.m",3,], bp.temp,
       angle = 90, length=.05, code = 3, lwd = 1)
bp.temp <- barplot(pars.arr["prop.exc.f",2,], xlim = c(0,1), names.arg = rep("",length(ds.nm)), horiz=2, border = NA,
                   las = 2, bty = "n", ylab = "", main="", col = cols, cex.axis = cex.names)
mtext("D", side = 3, line = .2, cex = cex.names)
arrows(pars.arr["prop.exc.f",1,], bp.temp,
       pars.arr["prop.exc.f",3,], bp.temp,
       angle = 90, length=.05, code = 3, lwd = 1)
## mtext("males", line = .5, side = 3, outer = T, cex = cex.names, at = .2)
## mtext("females", line = .5, side = 3, outer = T, cex = cex.names, at = .7)
mtext("proportional contribution of transmission by each route", side = 1, line = 0.2, outer = T, cex = cex.names, at = .4)
## legend("bottomright", c("before couple formation", "extra-couple",
##                 "within-couple"), ncol = 2, border = NA, fill = cols, cex = 1, bty = "n")
######################################################################
dev.off()

######################################################################
## Figure 3 Old: proportional contibution of each route for OBSERVED Couples
###################################################################### 
## plot
beside <- F
if(beside)
  {
    ords <- 3:1
    space <- c(.2,1)
  }else{
    ords <- 1:3
    space <- 0.1
  }
len <- .01
cex.labs <- .7
pdf("Figure S-Old 3-all observed breakdown.pdf", width = 4.21, height = 4.5)
par(oma = c(3.5,5.5,0,.5))
layout(matrix(c(1,2,4,6,1,3,5,7),nr=4,nc=2), widths = c(1,1), heights = c(.25,.8,.8,.8))
par(mar=rep(0,4))
plot(0,0, axes=F, type = "n")
cols <- c("dark gray","red", "dodgerblue1")[ords] # rev because reversing color order
legend("bottomright", c("before couple formation", "extra-couple",
                        "within-couple"), ncol = 2, border = NA, fill = cols, cex = 1, bty = "n")
par(mar=c(.5,1.5,.5,.5))
xlim <- c(0,1)
## discordant
bp.temp <- barplot(mbdn.dc3[ords,2,], beside = beside, names.arg = ds.nm.wa, las =2, space = space,
                   horiz = T, col = cols, xlim = xlim,  border = NA, axes = F)
mtext("males", side = 3, line = 0, cex = cex.labs)
mtext("serodiscordant", side = 2, line = 5.5, cex = cex.labs)
axis(1, seq(0,1,by=.2), labels = NA, las = 2)
if(beside) arrows(mbdn.dc3[ords,1,], bp.temp, mbdn.dc3[ords,3,], bp.temp,
                  angle = 90, length=len, code = 3, lwd = 1)
bp.temp <- barplot(fbdn.dc3[ords,2,], beside = beside,
                   names.arg = rep("", length(ds.nm.wa)), las =2, space = space,
                   horiz = T, col = cols, xlim = xlim,  border = NA, axes = F)
axis(1, seq(0,1,by=.2), labels = NA, las = 2)
if(beside) arrows(fbdn.dc3[ords,1,], bp.temp, fbdn.dc3[ords,3,], bp.temp,
                  angle = 90, length=len, code = 3, lwd = 1)
mtext("females", side = 3, line = 0, cex = cex.labs)
## concordant
bp.temp <- barplot(mbdn.cc[ords,2,], beside = beside, names.arg = ds.nm.wa, las =2, space = space,
                   horiz = T, col = cols, xlim = xlim,  border = NA, axes = F)
axis(1, seq(0,1,by=.2), labels = NA, las = 2)
mtext("concordant positive", side = 2, line = 5.5, cex = cex.labs)
if(beside) arrows(mbdn.cc[ords,1,], bp.temp, mbdn.cc[ords,3,], bp.temp,
                  angle = 90, length=len, code = 3, lwd = 1)
## mtext("C", side = 3, line = -.5, cex = .8)
bp.temp <- barplot(fbdn.cc[ords,2,], beside = beside,
                   names.arg = rep("", length(ds.nm.wa)), las =2, space = space,
                   horiz = T, col = cols, xlim = xlim,  border = NA, axes = F)
axis(1, seq(0,1,by=.2), labels = NA, las = 2)
if(beside) arrows(fbdn.cc[ords,1,], bp.temp, fbdn.cc[ords,3,], bp.temp,
                  angle = 90, length=len, code = 3, lwd = 1)
## mtext("D", side = 3, line = -.5, cex = .8)
## all
bp.temp <- barplot(mbdn.all[ords,2,], beside = beside, names.arg = ds.nm.wa, las =2, space = space,
                   horiz = T, col = cols, xlim = xlim,  border = NA, axes = F)
mtext("any couple", side = 2, line = 5.5, cex = cex.labs)
axis(1, seq(0,1,by=.2), las = 2)
if(beside) arrows(mbdn.all[ords,1,], bp.temp, mbdn.all[ords,3,], bp.temp,
                  angle = 90, length=len, code = 3, lwd = 1)
## mtext("E", side = 3, line = -.5, cex = .8)
bp.temp <- barplot(fbdn.all[ords,2,], beside = beside,
                   names.arg = rep("", length(ds.nm.wa)), las =2, space = space,
                   horiz = T, col = cols, xlim = xlim,  border = NA, axes = F)
axis(1, seq(0,1,by=.2), las = 2)
if(beside) arrows(fbdn.all[ords,1,], bp.temp, fbdn.all[ords,3,], bp.temp,
                  angle = 90, length=len, code = 3, lwd = 1)
## mtext("F", side = 3, line = -.5, cex = .8)
mtext("proportional contribution of transmission by each route", side = 1, line = 2.2, outer = T, cex = .8, at = .4)
######################################################################
dev.off()


######################################################################
## Model fit table
load(file = "~/Dropbox/disc mod backup/Couples Model Revision 120811/Output/p.val.vec.Rdata")
p.vals.d <- signif(p.val.vec,3)
## load(file = "~/Documents/R files/discordant couples/cloud final/no deaths p.val.vec.Rdata")
## p.vals.nd <- signif(p.val.vec,3)
## pframe <- data.frame(p.vals.d, p.vals.nd)
pframe <- data.frame(p.vals.d)#, p.vals.nd)
names(pframe) <- c("mortality")
write.csv(pframe, file = "pvals.csv")


if(do.extra )
  {
######################################################################
    ## Figure 1 area trajectory diagrams
######################################################################
    pdf("fig 1 traj.pdf",
        width = 3.5, height = 3.5)
    ## Load Zambian workspace to get those couples
    load("~/Documents/R files/discordant couples/sims/runs/Zambia-Deaths-20120824-09:31-rho0.5/files/workspace.Rdata")
    source("~/Dropbox/disc mod backup/Couples Model Revision 120811/R scripts & Data Files/pcalc5.R")
    col.vec <- colors()[c(89,438,544,395)]
    par(mfrow=c(2,2), mar = c(1,3,0,0), oma = c(0,0,0,0))
    last.year <- 2012
    cpls <- c(140:180,which(dat$ser<4)[c(50:90)]) # zam
    ## cpls <- c(1901,which(dat$ser<4)[c(38,5)]) # zim
    ## cpls <- c(147,388)
    cpls <- c(284, 233, 283, 501)
    for(cc in 1:length(cpls))
      {
        cpl <- cpls[cc]
        xticks <- T
        if(cc %in% 1:2) xticks <-  F
        ctraj.area(medpars, dat, browse =F, width = 4, height = 3, last.year = last.year,
                   plot.cpls = cpl, x.ticks = xticks, ylab = "", 
                   surv = survive,                   # show surv curves o plot
                   cols = col.vec, cex.t = .8,
                   nsurv = F, do.pdf = F,            # don't show marginal curves
                   dead = F, cprob = T,
                   survive = survive,                # plot point at survival curve
                   blob.cex = 1, show.age = F, show.susc = T, cex.leg = .8, do.leg = F, drp = 2,
                   t.lwd = 1.5,
                   pdf.name = "pos traj plots.pdf")
                                        #    mtext(cpl, side = 3, line = -3)
                                        #axis(1, seq(xlim[1],xlim[2], by = 60), seq(xlim[1],xlim[2], by = 60)/12 + 1900, pos = 0, las = 2,outer=F)
      }
    xlim <- c(80*12, (last.year-1900)*12)
                                        #mtext("probability in couple serostatus category / population prevalence", side = 2, line = -1, outer = T)
    dev.off()
###################################################################### 

    pdf("fig 1 leg.pdf", w = 3, h = 3)
    par(mar = rep(0,4))
    plot(0,0)
    legend("top", c("M-F-","M+F-","M-F+","M+F+"), fill =  col.vec, border = "black", bty = "n")
    dev.off()

    ## ######################################################################
    ## ## Figure 1 diagrams
    ## ######################################################################
    ## pdf("fig 1 traj.pdf",
    ##     width = 6.5, height = 4)
    ## ## Load Zambian workspace to get those couples
    ## load("~/Documents/R files/discordant couples/cloud final/deaths/Zambia-Deaths-20120608-12:05/files/workspace.Rdata")
    ## source("~/Documents/R files/discordant couples/reldurmod/pcalc5.R")
    ## par(mfrow=c(1,2), mar = c(4,5,0,0), oma = c(4,0,0,0))
    ## last.year <- 2012
    ## cpls <- c(140:180,which(dat$ser<4)[c(50:90)]) # zam
    ## ## cpls <- c(1901,which(dat$ser<4)[c(38,5)]) # zim
    ## cpls <- c(147,388)
    ## for(cc in 1:length(cpls))
    ##   {
    ##     cpl <- cpls[cc]
    ##     ctraj(medpars, dat, browse =F, width = 4, height = 3, last.year = last.year,
    ##           plot.cpls = cpl, x.ticks = F, ylab = "", 
    ##           surv = survive,                   # show surv curves o plot
    ##           nsurv = F, do.pdf = F,            # don't show marginal curves
    ##           dead = T, cprob = T,
    ##           survive = survive,                # plot point at survival curve
    ##           blob.cex = 2, show.age = F, show.susc = T, cex.leg = 1, do.leg = F, drp = 2,
    ##           t.lwd = 1.5,
    ##           pdf.name = "pos traj plots.pdf")
    ##     mtext(cpl, side = 3, line = -3)
    ##     axis(1, seq(xlim[1],xlim[2], by = 60), seq(xlim[1],xlim[2], by = 60)/12 + 1900, pos = 0, las = 2,outer=F)
    ##   }
    ## xlim <- c(80*12, (last.year-1900)*12)
    ## mtext("probability in couple serostatus category / population prevalence", side = 2, line = -1, outer = T)
    ## dev.off()
    ## ###################################################################### 

    ## ###################################################################### 
    ## ## Figure 99 - probability of a couple surviving to sampling given serostatus and age and rel history
    ## ######################################################################
    ## load("~/Documents/R files/discordant couples/cloud final/deaths/Zambia-Deaths-20120608-12:05/files/workspace.Rdata")
    ## source("~/Documents/R files/discordant couples/reldurmod/pcalc5.R")
    ## hh.log <- dat$ser == 1
    ## mm.log <- dat$ser == 2
    ## ff.log <- dat$ser == 3
    ## cex <- .8
    ## medpars <- pars[parnames,2]
    ## allstates <- pcalc(medpars, dat = dat, trace = T, give.pis=T, survive=survive, lrho.sd = lrho.sd)$allstates
    ## names(allstates)
    ## p.surv <- data.frame(mind = rep(NA, nrow(dat)), find = NA)
    ## p.surv[mm.log | hh.log,"mind"] <- rowSums(allstates[mm.log | hh.log,c("mb.A","me.A","hb1b2A","hbeA","hbpA","he1e2A")]) / rowSums(allstates[mm.log | hh.log,c("mb.","me.","hb1b2","hbe","hbp","he1e2")])
    ## p.surv[ff.log | hh.log,"find"] <- rowSums(allstates[ff.log | hh.log,c("f.bA","f.eA","hb2b1A","hebA","hpbA","he2e1A")]) / rowSums(allstates[ff.log | hh.log,c("f.b","f.e","hb2b1","heb","hpb","he2e1")])
    ## head(p.surv[mm.log,])
    ## head(p.surv[ff.log,])
    ## head(p.surv[hh.log,])

    ## par(mfrow=c(2,1), mar = rep(1,4))
    ## smoothScatter(dat$mage[mm.log|hh.log]/12, p.surv$mind[mm.log|hh.log],
    ##               cex = 1, col = "black",las = 2, cex.axis = 1,
    ##               main = "", xlab = "", ylab = "", nrpoints = 0,
    ##               xlim = c(0,60), ylim = c(0,1))
    ## smoothScatter(dat$fage[ff.log|hh.log]/12, p.surv$find[ff.log|hh.log],
    ##               cex = 1, col = "black",las = 2, cex.axis = 1,
    ##               main = "", xlab = "", ylab = "", nrpoints = 0,
    ##               xlim = c(0,60), ylim = c(0,1))

    ## ######################################################################
    ## pdf("Fig 6a - male index age vs surv.pdf", w = 3, h = 3)
    ## cols <- rep(NA, nrow(dat))
    ## cols[hh.log] <- "black"
    ## cols[mm.log] <- "dark green"
    ## cols[ff.log] <- "purple"
    ## scatterhist(dat$mage[mm.log|hh.log]/12, p.surv$mind[mm.log|hh.log], xlim = c(10,60), col = cols[mm.log|hh.log])
    ## dev.off()
    ## pdf("Fig 6b - female index age vs surv.pdf", w = 3, h = 3)
    ## scatterhist(dat$fage[ff.log|hh.log]/12, p.surv$find[ff.log|hh.log], xlim = c(10,60), col = cols[ff.log|hh.log])
    ## legend("topright", c("M+F-","M-F+","M+F+"), pch = 19, bty = "n", col = c("dark green", "purple", "black"))
    ## dev.off()




######################################################################
    ## Comparison of deaths to no deaths
    ## MUST HAVE RUN THIS SCRIPT FOR survive= T & = F already
######################################################################
######################################################################
    load(file = "~/Documents/R files/discordant couples/cloud final/pars.arr.Rdata")
    pars.d <- pars.arr
    ## Load no deaths pars.arr.Rdata (from run of all fits without mortality)
    load(file = "~/Documents/R files/discordant couples/cloud final/no deaths pars.arr.Rdata")
    pars.nd <- pars.arr

    rr.surv <- t(signif(pars.d[hazs,2,] / pars.nd[hazs,2,],2))
    rownames(rr.surv) <- ds.nm
    apply(pars.d[hazs[1:2],2,] / pars.nd[hazs[1:2],2,],2,range) #bef
    apply(pars.d[hazs[3:4],2,] / pars.nd[hazs[3:4],2,],2,range) #exc
    apply(pars.d[hazs[5:6],2,] / pars.nd[hazs[5:6],2,],2,range) #part

    write.csv(rr.surv[nrow(rr.surv):1,], file = "~/Documents/R files/discordant couples/cloud final/rr surv tab.csv")


######################################################################
    ## Show predicted proportion of transmission over the next year due to
    ## extracouple transmission for analyses ignoring death, used to argue
    ## the point that even though we ignore ARV effects on reduced
    ## mortality it's not a problem.
    load("~/Documents/R files/discordant couples/cf old 2/no deaths output/pars.arr.Rdata")
    pars.arr["prop.exc.m",,]
    pars.arr["prop.exc.f",,]
    prtab <- t(apply(pars.arr[c("prop.exc.m","prop.exc.f"),,],3,cis))
    rownames(prtab) <- rev(ds.nm.wa)
    prtab <- prtab[nrow(futtab):1,]
    write.csv(prtab, file="proportion incidence future breakdown No Deaths.csv")
  }
