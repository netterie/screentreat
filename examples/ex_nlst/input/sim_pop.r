#############################################################################
## File Name: sim_pop.r

## File Purpose: Simulate a study population based on
##    summary statistics reported in NLST trial
## Author: Leslie Mallinger
## Date: 3/20/13
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

n <- 2000
study.year <- 2003 
agetab <- data.frame(agemin=c(50, 60, 65, 70),
                     agemax=c(59, 64, 69, 74),
                     n=c(428, 306, 178, 88))
sextab <- data.frame(group=c(0, 1),
                     n=c(410, 590))
stagetab <- data.frame(group=c('I', 'II', 'III', 'IV'),
                       n=c(311, 79, 249, 361))

savename <- 'input_data.csv'

############################################################
# establish functions 
############################################################
addvar <- function(dset, varname, vartab) {
    d <- dset[sample(dset$id, nrow(dset), replace=FALSE), ]
    d[, varname] <- rep(vartab$group, vartab$n)
    d <- d[order(d$id), ]
    return(d)
}

############################################################
# simulate population 
############################################################
set.seed(seed)
pop <- data.frame(id=1:n,
                  agemin=rep(agetab$agemin, agetab$n),
                  agemax=rep(agetab$agemax, agetab$n))

pop$agegroup <- paste0(pop$agemin, '-', pop$agemax)
pop <- ddply(pop,
             .(agegroup),
             function(d) {
                 agemin <- unique(d$agemin)
                 agemax <- unique(d$agemax)
                 d$age <- round(runif(nrow(d), agemin, agemax))
                 return(d)
             })
pop <- subset(pop, select=-c(agemin, agemax, agegroup))
pop <- addvar(pop, 'male', sextab)
pop <- addvar(pop, 'stage', stagetab)

############################################################
# save data to file
############################################################
write.csv(pop, file=savename, row.names=FALSE)

