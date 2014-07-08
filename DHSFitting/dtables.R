library(data.table)

data.table()

n <- 10^3
k <- 9
serostates <- matrix(0,n,k)
serostates <- as.data.table(serostates)
setnames(serostates, 1:k, c('s..', 'mb.a1', 'mb.a2', 'mb.', 'f.ba1', 'f.ba2', 'f.b', 'hb1b2', 'hb2b1'))
serostates[, `:=`(s.. = 1)]
serostates

pre.couple <- function(serostates, sexually.active) {
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

sexually.active <- rbinom(n, 1,.5)==1

p.m.bef <- .5
p.f.bef <- .8

pre.couple(serostates, sexually.active)

head(serostates,50)
