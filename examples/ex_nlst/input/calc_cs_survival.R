#############################################################################
## File Name: user_options.r

## File Purpose: 
## Author: Leslie Mallinger
## Date: 8/5/2013 
## Edited on: 

## Additional Comments: 
#############################################################################

############################################################
# set up R 
############################################################
rm(list=ls())
library(foreign)
library(reshape)

############################################################
# establish functions 
############################################################
# simplify variable names
simpvarname <- function(dset, vars) {
    for (var in vars) {
        names(dset)[grep(var, names(dset))] <- var
    }
    return(dset)
}

# prep data
prep.data <- function(infile, by.var) {
    dset <- read.csv(infile, stringsAsFactors=FALSE)
    dset <- simpvarname(dset, 
                        c('Sex', 'stage', 'Interval', 
                          'SE', 'Lower', 'Upper'))
    dset <- subset(dset,
                   select=-c(Page.type, N, SE, Lower, Upper))
    dset <- rename(dset, c(Cause.Spec='Survival'))
    dset$Interval <- as.numeric(sub(' mo', '', dset$Interval))/12
    dset$Survival <- as.numeric(sub('%', '', dset$Survival))/100
    dset <- transform(dset, lambda=-log(Survival)/Interval)
    dset <- ddply(dset,
                  c(by.var),
                  summarize,
                  lambda=mean(lambda))
}

############################################################
# run calculations 
############################################################
# by stage
stage <- prep.data(infile='survival_by_stage.txt', by.var='stage')

# by sex
sex <- prep.data(infile='survival_by_sex.txt', by.var='Sex')


