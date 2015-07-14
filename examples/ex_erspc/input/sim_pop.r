#############################################################################
## File Name: sim_pop.r

## File Purpose: Simulate a study population based on
##    summary statistics
## Author: Leslie Mallinger
## Date: 5/13/2013
## Edited on: 

## Additional Comments: 
#############################################################################

############################################################
# set up R
############################################################
rm(list=ls())
library(xtable)
library(reshape)

############################################################
# set parameters 
############################################################
seed <- 98103

male <- 1

n <- 1000

study.year <- 2000

agetab <- data.frame(group=c('55-64', '65-69', '70-74'),
                     rx=0,
                     n=c(.5*n, .3*n, .2*n))
stagetab <- data.frame(group=c('L', 'R', 'D'),
                       rx=0,
                       n=c(.3*n, .3*n, .4*n))
txtab <- data.frame(group=c('RP', 'RT', 'CM'),
                    rx=0,
                    n=c(.2*n, .4*n, .4*n))

savename <- 'input_data.csv'

############################################################
# establish functions 
############################################################
makedistvar <- function(dset, varname, vartab) {
    var <- round(rnorm(nrow(dset), mean=vartab$mean, sd=vartab$sd))
    while (sum(var<vartab$min | var>vartab$max)>0) {
        var[var<vartab$min | var>vartab$max] <-
            round(rnorm(sum(var<vartab$min | var>vartab$max),
                        mean=vartab$mean,
                        sd=vartab$sd))
    }
    return(var)
}

adddistvar <- function(dset, varname, vartab) {
    var0 <- makedistvar(subset(dset, rx==0),
                        varname,
                        subset(vartab, rx==0))
    var1 <- makedistvar(subset(dset, rx==1),
                        varname,
                        subset(vartab, rx==1))
    dset[, varname] <- c(var0, var1)
    return(dset)
}

addcutvar <- function(dset,
                      varname.orig,
                      varname.new,
                      varbreaks,
                      varlabels) {
    dset[, varname.new] <- cut(dset[, varname.orig],
                               breaks=varbreaks,
                               labels=varlabels)
    return(dset)
}

addgroupvar <- function(dset, varname, vartab) {
    d0 <- subset(dset, rx==0)
    d1 <- subset(dset, rx==1)
    d <- rbind(d0[sample(d0$id, nrow(d0), replace=FALSE), ],
               d1[sample(d1$id-nrow(d0), nrow(d1), replace=FALSE), ])
    d[, varname] <- rep(vartab$group, vartab$n)
    d <- d[order(d$id), ]
    return(d)
}

dist_range_unif = function# Expand a range according to a uniform distribution

(n,
    ### Number of values to return from a uniform 
    ### distribution defined by rng
 rng, 
    ### String containing a range of values, separated by
    ### splitchar, to be expanded into n individual values
 splitchar='-',
    ### String character separating two values in rng
 to.round.int=FALSE
    ### If results will be rounded to the nearest integer,
    ### set to TRUE. This will expand the allowable range 
    ### for simulation so that integers in the range are
    ### evenly distributed after rounding 
) {

    rng.min=as.numeric(unlist(strsplit(rng, splitchar))[1])
    rng.max=as.numeric(unlist(strsplit(rng, splitchar))[2])
    if (to.round.int) {
        return(runif(n, rng.min-0.49999, rng.max+0.49999))
    } else {
        return(runif(n, rng.min, rng.max))
    }

### A vector of length n
}

expand_rangevar = function# Take in a data frame in which one column is a range variable, then expand it to integers according to a uniform distribution

(dset,
    ### Data frame containing range variable 
 rangevar, 
    ### Name of range variable in dset
 newvar,
    ### Name of new variable to be created with expanded values
 splitchar='-'
    ### String character separating values in rangevar
) {
    names(dset)[names(dset)==rangevar] = 'rangevar'
    dset = ddply(dset,
                  .(rangevar),
                  function(x) {
                      x[, newvar] = round(dist_range_unif(nrow(x),
                          unique(as.character(x$rangevar)),
                          splitchar,
                          to.round.int=TRUE))
                      return(x)
                  })
    names(dset)[names(dset)=='rangevar'] = rangevar
    return(dset)
}


############################################################
# simulate population 
############################################################
set.seed(seed)
pop <- data.frame(id=1:n,
                  male=male,
                  rx=0,
                  study_year=study.year)
pop <- addgroupvar(pop, 'agegroup', agetab)
pop <- expand_rangevar(dset=pop,
                       rangevar='agegroup',
                       newvar='age',
                       splitchar='-')
pop <- addgroupvar(pop, 'stage', stagetab)
pop <- addgroupvar(pop, 'tx', txtab)

pop <- subset(pop, select=-c(rx, agegroup, study_year))

############################################################
# save data to file
############################################################
write.csv(pop, file=savename, row.names=FALSE)

