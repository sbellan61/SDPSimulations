####################################################################################################
## This script takes the AIDSinfo data base and interpolates epidemic
## curves by gender for HIV prevalence, infectious HIV prevalence, and
## ART coverage (% of people living with HIV on ART). Interpolation is
## done on a logistic scale & then back-transformed.
## 
setwd('/home1/02413/sbellan/DHSFitting/')    # set working directory
load('data files/allDHSAIS.Rdata')      # load DHS data
## source('aidsinfo.R')

## Load AIDSinfo Database
un <- read.csv('data files/AIDSinfo_2013_en.csv')
## please DOWNLOAD 'WPP2012_DB04_POPULATION_ANNUAL.CSV' FROM http://esa.un.org/wpp/ASCII-Data/DISK_NAVIGATION_ASCII.htm
pop <- read.csv('data files/WPP2012_DB04_POPULATION_ANNUAL.CSV') # load population size data, 

##  rename some countries
un$Area <- as.character(un$Area)        
un$Area[un$Area=="United Republic of Tanzania"] <- 'Tanzania'
un$Area[grepl('Ivoire',un$Area)] <- 'Cote dIvoire'
un$Area[un$Area=='Democratic Republic of the Congo'] <- 'DRC'
un$Area <- factor(un$Area) ## change back to factor
for(ii in 1:5) un[,ii] <- factor(un[,ii]) ## set as factors
levels(un[,1])
show <- c('Area','Time.Period','Data.Value') # variables to show

## gather HIV prevalence (modelled)
psel <- un$Indicator=='HIV Prevalence'       # select  HIV prevalence
tail(table(un$Subgroup[psel])[order(table(un$Subgroup[psel]))],10) # number of measurements by subgroup
prev.mod <- psel & un$Subgroup=='Adults (15-49) estimate modelled' # select adults subgroup
tail(table(un$Area[prev.mod])[order(table(un$Area[prev.mod]))],10)
ptab <- un[prev.mod,show] ##  create prevalence table
ptab <- ptab[order(ptab$Area, ptab$Time.Period),]
head(ptab)
tail(ptab)
## gather adult men & women living with HIV (to break down by gender)
plhiv <- un$Indicator=="People living with HIV"
## all sexes adult
aplhiv <- plhiv & un$Subgroup=='All sexes Adults (15+) estimate'
apltab <- un[aplhiv,show]
apltab <- apltab[order(apltab$Area,apltab$Time.Period),]
head(apltab)
tail(apltab)
## females
fplhiv <- plhiv & un$Subgroup=='Females Adults (15+) estimate'
fpltab <- un[fplhiv,show]
fpltab <- fpltab[order(fpltab$Area,fpltab$Time.Period),]
head(fpltab)
tail(fpltab)
## do these have the same area-time labels?
sum(apltab$Area!=fpltab$Area | apltab$Time.Period!=fpltab$Time.Period)
sum(ptab$Area!=fpltab$Area | ptab$Time.Period!=fpltab$Time.Period)
## Yes, so put them in one data frame
tab <- data.frame(country = ptab$Area, year = ptab$Time.Period, prev = ptab$Data.Value,
                  aplhiv = apltab$Data.Value,
                  mplhiv = apltab$Data.Value - fpltab$Data.Value, 
                  fplhiv = fpltab$Data.Value)
tab$country <- factor(as.character(tab$country))
head(tab)
levels(tab$country)

## gather ART Coverage, note that it was only started in 2005
asel <- un$Indicator=="People receiving antiretroviral therapy"
tail(table(un$Subgroup[asel])[order(table(un$Subgroup[asel]))],10)
asel <- asel & un$Subgroup=='Total'
head(un[asel,],10)
artab <- un[asel,show]
artab <- artab[order(artab$Area,artab$Time.Period),]
head(artab)
tail(artab)
## more countries with ART estimates than prevalence
length(unique(artab$Area))
nlevels(tab$country)
## select subset with prevalence
artab <- artab[which(paste(artab$Ar, artab$Time) %in% paste(tab$country, tab$year)),]
##  match country years between ART coverage and other table
mch <- match(paste(artab$Ar, artab$Time), paste(tab$country, tab$year))
tab$art <- NA ## add ART coverage to table
tab$art[mch] <- artab$Data.Value
tab$art[is.na(tab$art)] <- 0
## Error in Senegal 2010 ART coverage, interpolate 2009 & 2011
tab$art[tab$country=='Senegal' & tab$year=='2010'] <- round(.5*(tab$art[tab$country=='Senegal' & tab$year=='2009'] + tab$art[tab$country=='Senegal' & tab$year=='2011']))


## Any countries missing?
unique(dat$epic.nm)[!unique(dat$epic.nm) %in% levels(tab$country)]
##  only keep countries for which we have DHS data
tab <- tab[tab$country %in% unique(dat$epic.nm),]
tab$country <- factor(as.character(tab$country)) ## re-factor to reduce number of levels
nlevels(tab$country)
dim(tab)

## Get population sizes for men & women 15-49 to use as denominators for gender specific-prevalences
head(pop)
names(pop)[2] <- 'country'
## Relabel countries to make sure they match between population & AIDSinfo data bases
pop$country <- as.character(pop$country)
pop$country[grepl('Tanzania', pop$country)] <- 'Tanzania'
pop$country[grepl('Democratic Republic of the Congo', pop$country)] <- 'DRC'
pop$country[grepl("d'Ivoire", pop$country)] <- 'Cote dIvoire'
pop$country <- factor(pop$country)
levels(pop$country)[52] <- 'Cote dIvoire'
levels(tab$country)[!levels(tab$country) %in% levels(pop$country)]
apop <- pop[pop$country %in% levels(tab$country),]
apop$country <- factor(as.character(apop$country))
mean(pop$country %in% unique(tab$country))
unique(apop$country)
head(apop)
apop$Value <- apop$Value*1000           # in K's

## Get m/f populations
tab$madpop <- NA
tab$fadpop <- NA
for(cc in 1:nlevels(apop$country)) {
  for(yy in 1990:2012) {                # sum up adult population size (15+)
    count <- levels(apop$country)[cc]
    sel <- apop$country == count & apop$Time == yy
    temp <- apop[sel,]
    tab$madpop[tab$year==yy & tab$country == count] <- sum(temp[temp$Sex=='Male' & temp$AgeGrpStart >14,'Value'])
    tab$fadpop[tab$year==yy & tab$country == count] <- sum(temp[temp$Sex=='Female' & temp$AgeGrpStart >14,'Value'])
  }
}

####################################################################################################
## Calculate prevalence values
tab$prev.all <- tab$aplhiv / (tab$madpop + tab$fadpop) # 15 + prevalence
tab$mprev.all <- (tab$mplhiv/tab$madpop) ## male prevalence is ppl living w/ HIV 15+ / # pppl in age group (from UNDP ests)
tab$fprev.all <- (tab$fplhiv/tab$fadpop) ## female prevalence is ppl living w/ HIV 15+ / # pppl in age group (from UNDP ests)
tab$mart <- tab$art*(tab$mplhiv/tab$aplhiv) ##  men on ART equals people on ART times proportion of people living with HIV that are men
tab$fart <- tab$art*(tab$fplhiv/tab$aplhiv) # women
tab$art.cov <- tab$art/tab$aplhiv ##   ART coverage is number of people with ART divided by number of people living with HIV
tab$prev.inf <- tab$prev.all * (1-tab$art.cov) ## infectious prevalence equals prevalence of people infected and not on ART
tab$mprev.inf <- (tab$mplhiv - tab$mart)/tab$madpop # proportion of adult men living w/ HIV & not on ART
tab$fprev.inf <- (tab$fplhiv - tab$fart)/tab$fadpop # proportion of adult women living w/ HIV & not on ART
head(tab)
tail(tab,5)

## only show columns that we're using
spectr <- tab[,c("year","country",'prev.all',"mprev.all","fprev.all",'prev.inf',"mprev.inf","fprev.inf","art.cov")]
                                        #pls <- c('epic.all','epicm.all','epicf.all','epic','epicm','epicf','art.cov')
pls <- c('prev.all','mprev.all','fprev.all','prev.inf','mprev.inf','fprev.inf','art.cov')
logit <- function(x) log(x/(1-x))
ilogit <- function(x) exp(x) / (1 + exp(x))
x.seq <- 1:((2013-1900)*12) ##country month code from 1900 to 2013 
for(pp in pls)  assign(pp, data.frame(cmc = x.seq)) ## initialize all prevalence labels
pdf(file.path('data files','prevalence interpolation models.pdf'), w = 8, h = 7) ## initialize plot
for(ii in 1:length(unique(spectr$country)))    { ## for each country
  cc <- unique(spectr$country)[ii] ## get country name
  temp <- spectr[spectr$country==cc,] ## create temporary country specific data set
  ## add 0 prevalence in 1980 to stabilize      
  temp <- rbind(data.frame(year = 1980, country = cc,
                           prev.all = 10^-6, mprev.all = 10^-6, fprev.all = 10^-6,
                           prev.inf = 10^-6, mprev.inf = 10^-6, fprev.inf = 10^-6,
                           art.cov = 10^-8), temp)
  temp$art.cov[temp$art.cov==0] <- 10^-8 ## reset any zero ART coverages to small values
  temp$cmc <- (temp$year -1900)*12
  par(mfrow = c(3,3), oma = c(0,0,1,0))
  for(pl in pls) { ## for each prevalence type
    tpl <- temp[,pl]
    plot(temp[,'year'], tpl, pch = 19, main = pl, xlim = c(1980,2011), ylim = c(0, .35)) ## plot data points
    tpl.lgt <- logit(tpl)##  take logistic transformation
    tpl.mod <- smooth.spline(tpl.lgt ~ temp$cmc) ## fit smooth spline model
    tpl.pred <- predict(tpl.mod, x.seq) ## interpolate model for all months
    tpl.pred$y <- ilogit(tpl.pred$y) ## take inverse logistic transformation
    assign(pl, data.frame(get(pl), tpl.pred$y)) ## add country data to data frame as a new column
    lines(1900+ tpl.pred$x/12, tpl.pred$y, col = "red") ## show interpolation on figure
  }
  mtext(cc, side = 3, line = -1, outer = T, adj = .5, cex = 1.2)
}
dev.off()

for(pl in pls) {
  temp <- get(pl)
  names(temp)[-1] <- as.character(unique(spectr$country))
  assign(pl, temp)
}

epicm <- mprev.inf
epicf <- fprev.inf
epicm.all <- mprev.all
epicf.all <- fprev.all

## name countries
## save all files
save(epicm, epicf, epicm.all, epicf.all,
     prev.all, mprev.all, fprev.all,
     prev.inf, mprev.inf, fprev.inf,
     art.cov, file = 'data files/epic.Rdata')
save(epicm, epicf, epicm.all, epicf.all,
     prev.all, mprev.all, fprev.all,
     prev.inf, mprev.inf, fprev.inf,
     art.cov, file = '/home1/02413/sbellan/SDPSimulations/data files/epic.Rdata')
