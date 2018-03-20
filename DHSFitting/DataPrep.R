rm(list=ls())
library(plyr)
if(Sys.info()['sysname']=='Darwin') setwd('~/Documents/R Repos/SDPSimulations/DHSFitting/')
## setwd('/home1/02413/sbellan/DHSFitting/')    # set working directory
#load("data files/alldhs raw.Rdata")         # DHS data
## load("data files/African_couplescombine.RData") # DHS + AIS
load("data files/epic.Rdata")               ## UNAIDS
load("data files/province.Rdata")           ## Province data
ccodes <- read.csv("data files/dhs country codes.csv")[,1:2] ## country codes
dirnm <- file.path('data files','DHS explore')
if(!file.exists(dirnm))  dir.create(dirnm)
## source('DataPrep.R')


####################################################################################################
pdf(file.path(dirnm,"DHS-AIS cleaning summary.pdf")) # initialize pdf showing errors/inconsistencies
####################################################################################################
raw <- Answers
## Rename some variables
names(raw) <- sub("_", "", names(raw))
while(sum(grepl("_", names(raw))>0))     names(raw) <- sub("_",".", names(raw))
while(sum(grepl("marriage", names(raw), ignore.case = T))>0)     names(raw) <- sub("marriage","mar", names(raw), ignore.case = T)
while(sum(grepl("marital", names(raw), ignore.case = T))>0)     names(raw) <- sub("marital","mar", names(raw), ignore.case = T)
while(sum(grepl("intercourse", names(raw), ignore.case = T))>0)     names(raw) <- sub("intercourse","intc", names(raw), ignore.case = T)
while(sum(grepl("Years", names(raw), ignore.case = T))>0)     names(raw) <- sub("Years","yr", names(raw), ignore.case = T)
while(sum(grepl("number.of.unions", names(raw), ignore.case = T))>0)     names(raw) <- sub("number.of.unions","num.un", names(raw), ignore.case = T)

####################################################################################################
## Organize into surveys with survey name, country code, & year of survey
names(raw)[names(raw)=='Countrycode.and.phase'] <- 'scode'
xtabs(~scode,raw)
raw$ccode <- substr(raw$scode,1,2)
raw$epic.nm <- ccodes[match(raw$ccode,ccodes[,1]),2]
## 
surveys <- unique(raw[,c('survey','scode','epic.nm')])
surveys$start <- NA
surveys$end <- NA
surveys$sname <- NA
for(ss in 1:nrow(surveys)) {
  surveys$start[ss] <- floor(min(raw$Minterview.cmc[raw$survey == surveys$survey[ss]]) /12 + 1900)
  surveys$end[ss] <- floor(max(raw$Minterview.cmc[raw$survey == surveys$survey[ss]])/12 + 1900)
  if(surveys$start[ss]==surveys$end[ss]) surveys$sname[ss] <- paste(surveys$epic.nm[ss], surveys$start[ss])
  else surveys$sname[ss] <- paste0(surveys$epic.nm[ss], ' ', surveys$start[ss],'-',substr(surveys$end[ss],3,4))
}
surveys <- surveys[,c('survey','scode','sname','epic.nm')]
surveys$sname[surveys$scode=='ET4'] <- 'Ethiopia 2005'
surveys$sname[surveys$scode=='ET6'] <- 'Ethiopia 2011'
surveys$epic.nm <- as.character(surveys$epic.nm)
raw$ds <- surveys$sname[match(raw$survey, surveys$survey)]
surveys

####################################################################################################
## Look at regions
ddply(raw, .(epic.nm), summarise, length(unique(Fregion)))
ddply(raw, .(epic.nm), summarise, length(unique(Mregion)))
ddply(raw, .(epic.nm), summarise, unique(Mregion))

ddply(raw, .(epic.nm), summarise, sum(is.na(Mregion)))
ddply(raw, .(epic.nm), summarise, sum(is.na(Fregion)))

## unique(rawn$survey)[!unique(rawn$survey) %in% unique(raw$survey)]
    
####################################################################################################
## Remove Malawi 2000, no HIV data
mal2000 <- raw$ds=='Malawi 2000'
raw <- raw[!mal2000,]
raw$ds <- factor(as.character(raw$ds))

####################################################################################################
## Summarize
xtabs(~ds,raw)
names(raw)
dim(raw)
raw <- data.frame(uid = 1:nrow(raw),raw)

####################################################################################################
## Match to epidemic curves
raw$epic.ind <- NA
raw$epic.nm <- NA
for(ii in 2:ncol(epicm)) {
    cc <- colnames(epicm)[ii]
    ind <- grepl(cc, raw$ds)
    if(sum(ind)>0) {
        raw$epic.ind[ind] <- ii
        raw$epic.nm[ind] <- cc
      }
  }
unique(raw[,c("ds","epic.ind","epic.nm")])
raw$group <- raw$epic.nm
w.africa <- c("Burkina Faso", "Cameroon", 'Cote dIvoire', "Ghana", "Guinea", "Liberia", "Mali", "Niger", "Senegal",
              "Sierra Leone")
raw$group <- raw$epic.nm
raw$group[raw$epic.nm %in% w.africa] <- "WA"
xtabs(~group, raw)
raw$group <- factor(raw$group)
nlevels(raw$group)

####################################################################################################
## Examine HIV test results: excluding 'noser'
xtabs(~FHIVresult + ds, raw)            # Senegal has a large proportion HIV2
xtabs(~MHIVresult + ds, raw)
xtabs(~FHIVresult + group, raw)  # WA has greatest proportion HIV2, but still rather small. Ignore strain differences.
xtabs(~MHIVresult + group, raw)
nas <-  c("ERROR : V-, W+, M+", "ERROR : V-, W+, M-", "ERROR : V-, W-, M+", # These are removed as if missing
          "Indeterminant", "Not enough samples to complete protocol")
noser <- is.na(raw$FHIVresult) | is.na(raw$MHIVresult)
noser <- noser | raw$FHIVresult %in% nas | raw$MHIVresult %in% nas
aggregate(noser, list(raw$group), function(x) signif(mean(x),3)*100) ## 

## create serostatuses
fp <- grepl("positive", raw$FHIVresult)
mp <- grepl("positive", raw$MHIVresult)
ser <- rep(NA, nrow(raw))
ser[fp & mp] <- 1 #"M+ F+"
ser[!fp & mp] <- 2 #"M+ F-"
ser[fp & !mp] <- 3 #"M- F+"
ser[!fp & !mp] <- 4 #"M- F-"
raw$ser <- ser
xtabs(~group + ser, raw)

####################################################################################################
## Look at clusters
par(mfrow=c(7,8), mar = c(2,1,2,.5), oma = c(2,0,0,0))
for(ii in 1:nlevels(raw$ds)) {
  hist(as.numeric(raw$MclusterId[raw$ds==levels(raw$ds)[ii]]), breaks = 0:3300,
                    xlab = '', main = levels(raw$ds)[ii], col = 'black')
  abline(v=60, col = 'red')
}
mtext('male cluster', side = 1, adj = .5, line = .5, outer = T, cex = 1)
par(mfrow=c(7,8), mar = c(2,1,2,.5), oma = c(2,0,0,0))
for(ii in 1:nlevels(raw$ds)) {
  hist(as.numeric(raw$FclusterId[raw$ds==levels(raw$ds)[ii]]), breaks = 0:3300,
                    xlab = '', main = levels(raw$ds)[ii], col = 'black')
  abline(v=60, col = 'red')
}
mtext('female cluster', side = 1, adj = .5, line = .5, outer = T, cex = 1)


####################################################################################################
## Exclude polygamous couples: excluding 'polyg'. Relying on male's accounting of it.
## 
## numWives is number of *other* wives for pre-2011 surveys & total # of wives for 2011 or later (we
## think based on the below)
xtabs(~ds + MnumWives, raw)[,1:5]             # how many have more than 1 current wife?
xtabs(~ds + FnumOtherWives, raw)[,1:5]             # how many have more than 1 current wife?
xtabs(~ds + is.na(MnumWives), raw)             # how many have more than 1 current wife?
mnw <- raw$MnumWives ## some surveys ask how many other wives, some ask how many total. Correct in this variable
to.change <- rownames(xtabs(~ds + MnumWives, raw))[xtabs(~ds + MnumWives, raw)[,1] ==0]
mnw[raw$ds %in% to.change] <- mnw[raw$ds %in% to.change] - 1
xtabs(~ds + mnw, raw)[,1:5]
cbind(xtabs(~group + mnw, raw)[,1:5], xtabs(~group + FnumOtherWives, raw)[,1:5])
polyg <- mnw > 0 | is.na(mnw) ## > 0 other wives?
xtabs(~ds + polyg, raw)
xtabs(~group + polyg, raw)
aggregate(polyg, list(raw$ds), function(x) signif(mean(x))*100)

####################################################################################################
## How does MnumWives compare to 'MNumber.of.wives.partners'?
xtabs(~MNumber.of.wives.partners + ds, raw)
xtabs(~is.na(MNumber.of.wives.partners) + ds, raw)
## nna <- !is.na(raw$MNumber.of.wives.partners) & !is.na(raw$MnumWives)
##plot(jitter(raw$MNumber.of.wives.partners, a = .1), jitter(raw$MnumWives, a = .1), xlim = c(0,5), ylim = c(0,5), main='comparison of # of wives variables')


####################################################################################################
## Exclude couples that aren't really couples: excluding 'ncpl'. Everyone should be since that's how
## they were merged.
xtabs(~ds + is.na(MmarStatus), raw)
xtabs(~ds + is.na(FmarStatus), raw)
xtabs(~ds + MmarStatus,raw)
xtabs(~ds + FmarStatus,raw)
mmar <- raw$MmarStatus %in% c("Married", "Living together", "Living with partner",NA) ## including NA since these are couples by definition!
fmar <- raw$FmarStatus %in% c("Married", "Living together", "Living with partner",NA)
mar <- mmar & fmar
aggregate(mar, list(raw$ds), function(x) signif(mean(x))*100)
ncpl <- !mar

####################################################################################################
## Exclude couples for which neither partner are in their first stable cohabiting partnership: Excluding 'no1scp'
## 
## Create variable indicating both are in their first & only marriage first need to make NA's
## unknowns to allow logical manipulation
levels(raw$Fnum.un) <- c(levels(raw$Fnum.un), "UNK")
levels(raw$Mnum.un) <- c(levels(raw$Mnum.un), "UNK")
raw$Mnum.un[is.na(raw$Mnum.un)] <- "UNK"
raw$Fnum.un[is.na(raw$Fnum.un)] <- "UNK"
## raw$Mnum.un[raw$ds=='Ethiopia 2005'] <- 'Once' # because they're all NAs? Nah, don't do it. Could screw up marital durations!
xtabs(~ds + Mnum.un, raw)
xtabs(~ds + Fnum.un, raw)
xtabs(~Fnum.un + Mnum.un, raw)
first.mar <- raw$Fnum.un=="Once" | raw$Mnum.un=="Once" 
no1scp <- !first.mar
xtabs(~ds + first.mar, raw)
aggregate(first.mar, list(raw$ds), function(x) signif(mean(x))*100)

####################################################################################################
## Exclude if interview of m & f was <= 1 month apart: 'int.gr1mon'
tint.diff <- (raw$Minterview.cmc - raw$Finterview.cmc)
int.gr1mon <- abs(tint.diff) > 1 | is.na(tint.diff)
aggregate(int.gr1mon, list(raw$ds), function(x) signif(mean(x),3)*100)

####################################################################################################
## Rename interview variable 'tint'
raw$tint <- raw$Minterview.cmc

####################################################################################################
## Calculaate ages using cmc (country-month-code; months since 1900)
mage <- (raw$Minterview.cmc-raw$Mdob.cmc)/12
fage <- (raw$Finterview.cmc-raw$Fdob.cmc)/12
## Compare with their variable
## par(mfrow=c(2,1))
## plot(raw$FcurrAge, fage)
## plot(raw$McurrAge, mage)
par(mfrow=c(4,4), mar = c(2,2,3,.5), oma = c(2,0,0,0))
for(ii in 1:nlevels(raw$group)) {
  hist(mage[raw$group==levels(raw$group)[ii]], breaks = 0:80,
                    xlab = '', main = levels(raw$group)[ii], col = 'black')
  abline(v=60, col = 'red')
}
mtext('male age', side = 1, adj = .5, line = .5, outer = T, cex = 1)
for(ii in 1:nlevels(raw$group)) {
  hist(fage[raw$group==levels(raw$group)[ii]], breaks = 0:80,
                    xlab = '', main = levels(raw$group)[ii], col = 'black')
  abline(v=50, col = 'red')
}
mtext('female age', side = 1, adj = .5, line = .5, outer = T, cex = 1)
too.old <- mage > 59 | fage > 49 ## these are the max ages for DHS data
aggregate(too.old, list(raw$ds), function(x) signif(mean(x),3)*100)        ## looks like Mozambique 2009 & Uganda 2011 didn't follow protocol

####################################################################################################
## Examine
###################################################################### 
## Calculate relationship durations & check that when both partners are in their first marriage,
## that they report consistent marital durations.
mmardur <- (raw$Minterview.cmc-raw$MDate.of.first.mar)/12
fmardur <- (raw$Finterview.cmc-raw$FDate.of.first.mar)/12
## Are these consistent with DHS' years since first marriage?
for(ii in 1:nlevels(raw$group)) hist( 12*(mmardur - raw$Myr.since.first.mar)[raw$group==levels(raw$group)[ii]], breaks = -1:13,
                    xlab = '', main = levels(raw$group)[ii], col = 'black')
mtext("our male marriage duration - DHS male marriage duration (months)", side = 1, adj = .5, line = .5, outer = T, cex = 1)
for(ii in 1:nlevels(raw$group)) hist( 12*(fmardur - raw$F.yr.since.first.mar)[raw$group==levels(raw$group)[ii]], breaks = -1:13,
                    xlab = '', main = levels(raw$group)[ii], col = 'black')
mtext("our female marriage duration - DHS female marriage duration (months)", side = 1, adj = .5, line = .5, outer = T, cex = 1)
## Yes, everything is < 12 months (which means they're just rounding to nearest year)

mardiff <- mmardur-fmardur## Are the male & female reports consistent with each other?
for(ii in 1:nlevels(raw$group)) hist(mardiff[raw$group==levels(raw$group)[ii]], breaks = -50:50,
                    xlab = '', main = levels(raw$group)[ii], col = 'black')
mtext("male marriage duration - female marriage duration (yrs)", side = 1, adj = .5, line = .5, outer = T, cex = 1)
## Create a single marital duration variable using the partner in their 1st partnership or average
## of both if both are in their 1st
ffirst <- raw$Fnum.un=="Once" & (!raw$Mnum.un=="Once" | is.na(raw$MDate.of.first.mar))
mfirst <- raw$Mnum.un=="Once" & (!raw$Fnum.un=="Once" | is.na(raw$FDate.of.first.mar))
bfirst <- raw$Mnum.un=="Once" & raw$Fnum.un=="Once" & !is.na(raw$MDate.of.first.mar) & !is.na(raw$FDate.of.first.mar)
raw$both.first.scp <- raw$Mnum.un=="Once" & raw$Fnum.un=="Once"
raw$mfun <- raw$Mnum.un=='Once'
raw$ffun <- raw$Fnum.un=='Once'
raw$mardur[ffirst] <- fmardur[ffirst]
raw$mardur[mfirst] <- mmardur[mfirst]
raw$mardur[bfirst] <- .5*(mmardur+fmardur)[bfirst]
raw$mardur.mon <- round(raw$mardur*12)  # must round months since averaging m & f month durations
raw$tmar <- raw$tint - raw$mardur.mon    # imputed cmc of marriage using mardur above

####################################################################################################
## % error in marriage duration for couples where both are in first marriage: Exclude 'mardur.25'
pmardiff <-  abs(mardiff/raw$mardur)
pmardiff[!bfirst] <- NA
range(pmardiff[bfirst],na.rm=T)
for(ii in 1:nlevels(raw$group)) {
  hist(100*pmardiff[bfirst & raw$group==levels(raw$group)[ii]], breaks = seq(0,200, by = 10),
                    xlab = '', main = levels(raw$group)[ii], col = 'black')
  abline(v=25, col = 'red')
}
mtext('% absolute difference in marriage duration listed by m & f', side = 1, adj = .5, line = .5, outer = T, cex = 1)
for(ii in 1:nlevels(raw$group)) {
  hist(mardiff[pmardiff > .25 & bfirst & raw$group==levels(raw$group)[ii]], breaks = -50:50,
                    xlab = '', main = levels(raw$group)[ii], col = 'black')
  abline(v=25, col = 'red')
}
mtext('mmardur-fmardur amongst couples where this was > 25% of mardur & both reported it', side = 1, outer = T, line = .5)
## exclude if difference is > 25% of average
mardur.25 <- pmardiff > .25
mardur.25[is.na(mardur.25)] <- F        # from 0 duration marriages
mardur.25 <- mardur.25 & bfirst         # only worried about error when both in first
mardur.25[is.na(mardur.25)] <- F
aggregate(mardur.25, list(raw$ds), function(x) signif(mean(x),3)*100)
## ## exclude if difference is > 25% of average & > 3 months: Doesn't improve things much!
## mardur.25a <- pmardiff > .25  & abs(mardiff) > (3/12)
## mardur.25a[is.na(mardur.25a)] <- F        # from 0 duration marriages
## mardur.25a <- mardur.25a & bfirst         # only worried about error when both in first
## mardur.25a[is.na(mardur.25a)] <- F
## xtabs(~ds + mardur.25a, raw)

for(ii in 1:nlevels(raw$group)) hist(1900 + raw$tmar[raw$group==levels(raw$group)[ii]]/12, breaks = seq(1960,2015, by = 1),
                    xlab = '', main = levels(raw$group)[ii], col = 'black')
mtext('couple formation date', side = 1, adj = .5, line = .5, outer = T, cex = 1)
xtabs(~group + is.na(mardur.mon), raw, subset = first.mar)
## No marital date?: exclude 'nomdate' (this exludes more than no1scp)
nomdate <- is.na(raw$mardur)
aggregate(nomdate, list(raw$ds), function(x) signif(mean(x),3)*100)


####################################################################################################
## Exclude if missing age at first intercourse of if flagged as
## inconsistent with FAge.at.first.intc.impFlag: 'no.aafi'
raw$MAge.at.first.intc.impFlag <- factor(raw$MAge.at.first.intc.impFlag)
raw$FAge.at.first.intc.impFlag <- factor(raw$FAge.at.first.intc.impFlag)
## First note that Uganda 2004-05 has all NA's for the imputed variable, so let's impute ourselves:
## If aafi > (aafmar + 1), set to NA
##
## What proportion are NAs first off?
sel <- raw$ds == 'Uganda 2004-05'
## F
mean(is.na(raw$FAge.at.first.intc[sel]))*100
raw$FAge.at.first.intc.imputed[sel] <- raw$FAge.at.first.intc[sel]
fsdam <- raw$FAge.at.first.intc > (raw$FAge.at.first.mar + 1) # f sexual debut after marriage (fsdam)
sum(fsdam[sel],na.rm=T)/sum(sel)*100 # 27% reported sexual debuts after marriage in this data set??
raw$FAge.at.first.intc.imputed[sel & fsdam] <- NA
## M
mean(is.na(raw$MAge.at.first.intc[sel]))*100
raw$MAge.at.first.intc.imputed[sel] <- raw$MAge.at.first.intc[sel]
msdam <- raw$MAge.at.first.intc > (raw$MAge.at.first.mar + 1) # m sexual debut after marriage (msdam)
sum(msdam[sel],na.rm=T)/sum(sel)*100 # 17% reported sexual debuts after marriage in this data set??
raw$MAge.at.first.intc.imputed[sel & raw$MAge.at.first.intc > (raw$MAge.at.first.mar + 1)] <- NA
## update flags to be No flag or NA if we flagged it
raw$MAge.at.first.intc.impFlag[sel] <- c("No flag", NA)[as.numeric(is.na(raw$MAge.at.first.intc.imputed[sel]))+1]
raw$FAge.at.first.intc.impFlag[sel] <- c("No flag", NA)[as.numeric(is.na(raw$FAge.at.first.intc.imputed[sel]))+1]
## Look at all countries
xtabs(~raw$FAge.at.first.intc.impFlag, raw)
xtabs(~raw$MAge.at.first.intc.impFlag, raw)
levels(raw$MAge.at.first.intc.impFlag) %in% levels(raw$FAge.at.first.intc.impFlag) # same levels for genders?
oks <- c("No flag", "After concep < 1 yr", "After conception < 1 year", "After marriage") # if after marriage by < 1 yr, make it at marriage (later in code)
fflag <- ! raw$FAge.at.first.intc.impFlag %in% oks ## This also deals with NA's
mflag <- ! raw$MAge.at.first.intc.impFlag %in% oks ## This also deals with NA's
aggregate(fflag, list(raw$ds), function(x) signif(mean(x),3)*100) # all Uganda 2004-5 is flagged?
aggregate(mflag, list(raw$ds), function(x) signif(mean(x),3)*100) # all Uganda 2004-5, & Malawi 2000 are flagged? not worrying about the latter since we don't have HIV data
## DHS' imputation removed nearly all of the >90 codes, so let's just
## remove the rest instead of trying to interpret them.
oflag <- raw$FAge.at.first.intc.imputed >90 | raw$MAge.at.first.intc.imputed >90
oflag[is.na(oflag)] <- T
no.aafi <- fflag | mflag | oflag
xtabs(~raw$MAge.at.first.intc.imputed, subset = !no.aafi)
xtabs(~raw$FAge.at.first.intc.imputed, subset = !no.aafi)
aggregate(no.aafi, list(raw$ds), function(x) signif(mean(x),3)*100)

####################################################################################################
## Calculate years of sexual activity before marriage/couple formation
######################################################################
## Who (among those left after DHS & my imputation of sexual debuts) reported sexual debuts after first marriage?
aggregate(raw$FAge.at.first.intc.imputed[!no.aafi]>(raw$FAge.at.first.mar)[!no.aafi], list(raw$ds[!no.aafi]), function(x) signif(mean(x,na.rm=T))*100)
aggregate(raw$MAge.at.first.intc.imputed[!no.aafi]>(raw$MAge.at.first.mar)[!no.aafi], list(raw$ds[!no.aafi]), function(x) signif(mean(x,na.rm=T))*100)
## More than 1 year after?
aggregate(raw$FAge.at.first.intc.imputed[!no.aafi]>(raw$FAge.at.first.mar+1)[!no.aafi], list(raw$ds[!no.aafi]), function(x) signif(mean(x,na.rm=T))*100)
aggregate(raw$MAge.at.first.intc.imputed[!no.aafi]>(raw$MAge.at.first.mar+1)[!no.aafi], list(raw$ds[!no.aafi]), function(x) signif(mean(x,na.rm=T))*100)
## Still gotta remove a chunk of data from these people!
ncode <- raw$MAge.at.first.intc.imputed < 90
for(ii in 1:nlevels(raw$group)) {
    hist((raw$MAge.at.first.mar - raw$MAge.at.first.intc.imputed)[ncode & raw$group==levels(raw$group)[ii]],
         main = levels(raw$group)[ii], xlab = "", col="black", breaks =-90:90)
    abline(v=0, col = "red")
}
ncode <- raw$FAge.at.first.intc.imputed < 90
mtext('(male) age at first marriage - age at first intercourse (yrs)', side = 1, adj = .5, line = .5, outer = T, cex = 1)
for(ii in 1:nlevels(raw$group)) {
    hist((raw$FAge.at.first.mar - raw$FAge.at.first.intc.imputed)[ncode & raw$group==levels(raw$group)[ii]],
         main = levels(raw$group)[ii], xlab = "", col="black", breaks =-90:90)
    abline(v=0, col = "red")
}
mtext('(female) age at first marriage - age at first intercourse (yrs)', side = 1, adj = .5, line = .5, outer = T, cex = 1)

#################################################################################################### 
## Calculate CMC of sexual debuts
raw$fysa <- fage - raw$FAge.at.first.intc.imputed # female years of sexual activity (fysa)
raw$mysa <- mage - raw$MAge.at.first.intc.imputed
raw$tfs <- raw$tint - raw$fysa*12       # time of female sexual debut (time female sex, tfs)
raw$tms <- raw$tint - raw$mysa*12       # male
## For individuals that reported age at sexual debut after the marriage date but by < 1 yr, set
## their sexual debut to their month of marriage
aggregate((raw$tms>(raw$tmar))[!no.aafi], list(raw$ds[!no.aafi]), function(x) signif(mean(x,na.rm=T))*100)
aggregate((raw$tfs>(raw$tmar))[!no.aafi], list(raw$ds[!no.aafi]), function(x) signif(mean(x,na.rm=T))*100)
aggregate((raw$tms>(raw$tmar+12))[!no.aafi], list(raw$ds[!no.aafi]), function(x) signif(mean(x,na.rm=T))*100)
aggregate((raw$tfs>(raw$tmar+12))[!no.aafi], list(raw$ds[!no.aafi]), function(x) signif(mean(x,na.rm=T))*100)
## Remove those >1 year apart
mar.bef.sex <- raw$tms > (raw$tmar + 12) | raw$tfs > (raw$tmar + 12)
## Set those within 1 year to marriage date
msel <- which(raw$tms>raw$tmar & !raw$tms>(raw$tmar+12))
raw$tms[msel] <- raw$tmar[msel]
fsel <- which(raw$tfs>raw$tmar & !raw$tfs>(raw$tmar+12))
raw$tfs[fsel] <- raw$tmar[fsel]
mar.bef.sex[is.na(mar.bef.sex)] <- F
aggregate(mar.bef.sex, list(raw$ds), function(x) signif(mean(x))*100)

graphics.off()


####################################################################################################
## people married too young?: 'mar.und.8'
for(ii in 1:nlevels(raw$group)) hist((fage-fmardur)[raw$group==levels(raw$group)[ii]], breaks = -5:100, xlab = "", main = levels(raw$group)[ii], col = 'black')
mtext('male age at first marriage', side = 1, adj = .5, line = .5, outer = T, cex = 1)
for(ii in 1:nlevels(raw$group)) hist((mage-mmardur)[raw$group==levels(raw$group)[ii]], breaks = -5:100, xlab = "", main = levels(raw$group)[ii], col = 'black')
mtext('female age at first marriage', side = 1, adj = .5, line = .5, outer = T, cex = 1)
## Data to remove: people that say they married under 10 ys old if they are in their first marriage
## age in months at interview
raw$fage <- (raw$Finterview.cmc - raw$Fdob.cmc)
raw$mage <- (raw$Minterview.cmc - raw$Mdob.cmc)
mar.und.8 <- raw$mardur.mon >= raw$mage-8*12 |raw$mardur.mon >= raw$fage-8*12
mar.und.8[is.na(mar.und.8)] <- F
head(raw[mar.und.8,c("ds","mage","fage","tint","tmar")])
aggregate(mar.und.8, list(raw$ds), function(x) signif(mean(x))*100)

## no one should be having sex < 5 yrs old: 'early.sex'
early.sex <- raw$tint-raw$tms +60>= raw$mage
early.sex <- early.sex | raw$tint-raw$tfs +60>= raw$fage
early.sex[is.na(early.sex)] <- F
raw[early.sex,c('ds',"tint","tms","mage","tfs","fage")]
aggregate(early.sex, list(raw$ds), function(x) signif(mean(x))*100)

####################################################################################################
## Correct for Ethiopia's different calendar
## 
## Ethiopia 2005: "The survey was fielded from April 27 to August 30,
## 2005." (p. 11, Ethiopia DHS 2005 report)
## difference in dates is then
diff2005 <- 105*12 + 4 - min(raw[raw$ds=="Ethiopia 2005","tint"])
raw[raw$ds=="Ethiopia 2005",c("tms","tfs","tmar","tint")] <- raw[raw$ds=="Ethiopia 2005",c("tms","tfs","tmar","tint")] + diff2005
## Ethiopia 2011:"Rawa collection took place over a five-month period
## from 27 December 2010 to 3 June 2011." (p. 10, Ethiopia DHS 2011
## report)
diff2011 <- 110*12 + 12 - min(raw[raw$ds=="Ethiopia 2011","tint"])
raw[raw$ds=="Ethiopia 2011",c("tms","tfs","tmar","tint")] <- raw[raw$ds=="Ethiopia 2011",c("tms","tfs","tmar","tint")] + diff2011

## Couples removed due to inconsistent data
errs <- mar.und.8 | mardur.25 | mar.bef.sex | early.sex |  int.gr1mon | too.old
## Couples removed due to missing data (other than HIV serostatus)
mis <- nomdate | no.aafi
## Missing HIV serostatus
aggregate(noser, list(raw$ds), function(x) signif(mean(x))*100)
## removed from analysis
rem <- errs | mis | noser | polyg
anal <- !rem
## removed only because they didn't have serostatus (for sensitivity analyses later comparing
## couples w/ & w/o serostatus
remser <- noser & !(mis | polyg | errs)
rem <- rem | grepl("Sao Tome", raw$ds) # don't have epidemic curve for sao tome
names(raw)[names(raw)=='Mregion'] <- 'region'
## Remove any countries with less than 100 couples left
show <- c("uid","ds","ser","tms","tfs","tmar","tint","mardur.mon","mage","fage", 'both.first.scp', # "circ",
          "epic.ind","epic.nm","group",'region') #, "m.k.arv","f.k.arv", "mlsp", "flsp", "mevtest","fevtest")
## extra variables for dissolution analysis
show.dis <- c("uid","ds","ser","tms","tfs","tmar","tint","mardur.mon","mage","fage", 'both.first.scp', # "circ",
          "epic.ind","epic.nm","group",'region','mfun','ffun') #, "m.k.arv","f.k.arv", "mlsp", "flsp", "mevtest","fevtest")
## show <- c("uid","ds","ser","tms","tfs","tmar","tint","mardur.mon","mage","fage", 'both.first.scp', # "circ",
##           "epic.ind","epic.nm","group") #, "m.k.arv","f.k.arv", "mlsp", "flsp", "mevtest","fevtest")
## ## extra variables for dissolution analysis
## show.dis <- c("uid","ds","ser","tms","tfs","tmar","tint","mardur.mon","mage","fage", 'both.first.scp', # "circ",
##           "epic.ind","epic.nm","group",'mfun','ffun') #, "m.k.arv","f.k.arv", "mlsp", "flsp", "mevtest","fevtest")
dat <- raw[!rem,show]
dat.dis <- raw[!rem,show.dis]
dim(dat)
## Check for any remaining NA's?
apply(dat,2,function(x) sum(is.na(x)))

####################################################################################################
## error table
crit <- cbind(noser, nomdate, no.aafi, polyg, too.old, int.gr1mon, mardur.25, mar.bef.sex, early.sex, mar.und.8)
nms <- c('survey','no HIV \ntest', 'no marital \nduration', 'no sexual \ndebut', 'polygamous', 'too old', 'interview \ndates differ \nby >1 month',
         'partner marital \ndurations differ \nby >25% of \n average', 'married >1yr \nbefore sex debut', 'sexual debut \n<5yrs old',
         'married under \n8yrs old')
remtab <- aggregate(crit, list(raw$ds), function(x) signif(mean(x),3)*100) ## 
names(remtab) <- nms
write.csv(remtab, file = file.path(dirnm, 'percentage disqualified by criteria (survey).csv'))
remtab.gr <- aggregate(crit, list(raw$group), function(x) signif(mean(x),3)*100) ## 
names(remtab.gr) <- nms
write.csv(remtab.gr, file = file.path(dirnm, 'percentage disqualified by criteria (group).csv'))


####################################################################################################
## error table: due to just this criteria
onlycrit <- crit
for(vv in 1:ncol(crit)) onlycrit[,vv] <- crit[,vv] & !apply(crit[,-vv], 1, any)
nms <- c('survey','no HIV \ntest', 'no marital \nduration', 'no sexual \ndebut', 'polygamous', 'too old', 'interview \ndates differ \nby >1 month',
         'partner marital \ndurations differ \nby >25% of \n average', 'married >1yr \nbefore sex debut', 'sexual debut \n<5yrs old',
         'married under \n8yrs old')
onlyremtab <- aggregate(onlycrit, list(raw$ds), function(x) signif(mean(x),3)*100) ## 
names(onlyremtab) <- nms
write.csv(onlyremtab, file = file.path(dirnm, 'percentage disqualified ONLY by this criteria (survey).csv'))
onlyremtab.gr <- aggregate(onlycrit, list(raw$group), function(x) signif(mean(x),3)*100) ## 
names(onlyremtab.gr) <- nms
write.csv(onlyremtab.gr, file = file.path(dirnm, 'percentage disqualified ONLY by this criteria (group).csv'))

####################################################################################################
## Start table recording couples excluded & why
ltab <- data.frame(n = xtabs(~group, raw), nohivtest = NA, polygamy = NA,
                   missing = NA, errors = NA, analyzed = NA)
ltab$nohivtest <- xtabs(~group + noser, raw)[,2]## update inclusion table
ltab$polygamy <- xtabs(~group + polyg, raw)[,2]
ltab$missing <- xtabs(~group + mis, raw)[,2]
ltab$errors <- xtabs(~group + errs, raw)[,2]
ltab$analyzed <- xtabs(~group + anal, raw)[,2]
## add # of couples where both are in 1st cohabitation in analyzed group
ltab$bothfirstSCP <- xtabs(~group + both.first.scp, dat)[,2]
serstat <- xtabs(~group + ser, dat)[,c(4,2,3,1)]
colnames(serstat) <- c('ccn','msd','fsd','ccp')
ltab <- data.frame(ltab,serstat)
serstat <- xtabs(~group + ser, raw)[,c(4,2,3,1)]
colnames(serstat) <- paste(c('ccn','msd','fsd','ccp'), 'raw')
ltab <- data.frame(ltab,serstat)

for(ii in 1:nrow(ltab)) {
    ltab[ii,8:12] <- paste(as.numeric(ltab[ii,8:12]),
                          " (", round(as.numeric(ltab[ii,8:12])/ltab[ii,7],3)*100, "%)", sep="")
  }
for(ii in 1:nrow(ltab)) {
    ltab[ii,3:7] <- paste(as.numeric(ltab[ii,3:7]),
                          " (", round(as.numeric(ltab[ii,3:7])/ltab[ii,2],3)*100, "%)", sep="")
  }
for(ii in 1:nrow(ltab)) {
    ltab[ii,13:16] <- paste(as.numeric(ltab[ii,13:16]),
                          " (", round(as.numeric(ltab[ii,13:16])/sum(as.numeric(ltab[ii,13:16])),3)*100, "%)", sep="")
  }
ltab
write.csv(ltab, file = file.path(dirnm, "loss table.csv"))
show <- c("uid","ds","ser","tms","tfs","tmar","tint","mardur.mon","mage","fage", # 'both.first.scp', # "circ",
          "epic.ind","epic.nm","group",'region') #, "m.k.arv","f.k.arv", "mlsp", "flsp", "mevtest","fevtest")
dat <- dat[,show]


## chi sq for difference in  serostatus b/w those that were removed vs not
allraw <- data.frame(raw, remser, noser, rem)
allraw$ser[noser] <- NA
save(allraw, file = file.path('data files', 'allRawDHSAIS.Rdata'))
save(allraw, file = file.path('/home1/02413/sbellan/SDPSimulations/data files/', 'allRawDHSAIS.Rdata'))
     
dat$uid <- 1:nrow(dat)
xtabs(~group, dat)
xtabs(~ds, dat)
nrow(dat)
dat$group <- factor(dat$group)
dat$ds <- factor(dat$ds)
save(dat, file = file.path('data files', 'allDHSAIS.Rdata'))
save(dat, file = file.path('/home1/02413/sbellan/SDPSimulations/data files/', 'allDHSAIS.Rdata'))

dat.dis$group <- factor(dat.dis$group)
dat.dis$ds <- factor(dat.dis$ds)
save(dat.dis, file = file.path('data files', 'allDHSAISdis.Rdata'))
save(dat.dis, file = file.path('/home1/02413/sbellan/SDPSimulations/data files/', 'allDHSAISdis.Rdata'))

ds.nm <- levels(dat$group)
save(ds.nm, file = file.path('data files', 'ds.nm.all.Rdata'))
save(ds.nm, file = file.path('/home1/02413/sbellan/SDPSimulations/data files/', 'ds.nm.all.Rdata'))

sum(dat$tms>dat$tmar)
sum(dat$tfs>dat$tmar)
sum(dat$tmar>dat$tint)
apply(dat,2,function(x) sum(is.na(x)))

survs <- data.frame(group = levels(dat$group), couples = as.numeric(xtabs(~group, dat)), surveys = NA)
for(cc in 1:nrow(survs)) survs$surveys[cc] <- paste(unique(dat$ds[dat$group==levels(dat$group)[cc]]), collapse = ', ')
write.csv(survs, file='/home1/02413/sbellan/DHSFitting/data files/DHS explore/samp sizes.csv')
