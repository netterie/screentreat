#############################################################################
## File Name: report_prep.r 

## File Purpose: 
## Author: Leslie Mallinger
## Date: 8/6/2013
## Edited on: 

## Additional Comments: 
#############################################################################

############################################################
# set up R 
############################################################
rm(list=ls())
setwd('~/projects/cantrance/screening/ex_nlst/reports')
library(plyr)
library(xtable)

# functions
str2df = function(str, ...) {
    if (str!="") {
        df = read.csv(con <- textConnection(str),
                      stringsAsFactors=FALSE,
                      ...)
        close(con)
    } else df = NA
    return(df)
}

str2vec = function(str, vectype='character', ...) {
    vec = scan(con <- textConnection(str),
               what=vectype,
               sep=',',
               ...)
    close(con)
    return(vec)
}

compile_list = function(str1, str2, str3) {
    l = list(str2df(str1),
             str2df(str2),
             str2df(str3))
    if (sum(is.na(l))==3) l <- NULL else l <- l[!is.na(l)]
    return(l)
}

get_lambda = function(param, values, k=NA) {
    return(switch(param,
                  median = log(2)/values,
                  mean = 1/values,
                  ksurv = -log(values)/k,
                  rate = values))
}

convert_HR = function(lst=NULL, surv_param, surv_k=NA, baseline_rate) {
    if (!is.null(lst)) {
        return(lapply(lst,
                      function(x) {
                          if ('stat'%in%names(x)) {
                              x$rate = get_lambda(param=surv_param,
                                                   values=x$stat,
                                                   k=surv_k)
                              x = transform(x, HR=rate/baseline_rate)
                              x = subset(x, select=-c(stat, rate))
                          } else x
                      }))
    } else return(lst)
}



############################################################
# prep population characteristics 
############################################################
# source code
source ('../input/user_options.R')

# times
times <- str2vec(times)

# helper data
rxname <- data.frame(rx=c(0, 1),
                     Mode=c('Clinical', 'Screened'))

# stage table
stgtab <- str2df(categorical_chars1)
stgtab <- rename(stgtab, c(stage='Clinical Stage',
                           prop='Proportion'))

# age table
agetab <- str2df(categorical_chars2)
agetab <- rename(agetab, c(agegroup='Age Group',
                           prop='Proportion'))

# sex table
sextab <- str2df(categorical_chars3)
sextab <- rename(sextab, c(male='Sex',
                           prop='Proportion'))
sextab$Sex <- revalue(as.character(sextab$Sex),
                      c('0'='Female',
                        '1'='Male'))


############################################################
# prep user options 
############################################################
# population stage distribution in presence of screening
scrtab <- str2df(scr_stg_dist)
scrtab <- rename(scrtab,
                 c(order='Order',
                   stage='Screened Stage',
                   prop='Proportion'))

# type of parameter for mortality
morttype <- switch(mort_param,
                  rate='cause-specific survival rate',
                  mean='mean cause-specific survival',
                  median='median cause-specific survival',
                  ksurv=paste0(mort_k,
                               '-year cause-specific survival'))
mort_rate <- get_lambda(param=mort_param,
                        values=mort_value,
                        k=mort_k)

# covariate-specific survival statistics for cause-specific mortality
morthrlist <- compile_list(mort_covar1, mort_covar2, mort_covar3)
morthrlist <- convert_HR(lst=morthrlist,
                         surv_param=mort_param,
                         surv_k=mort_k,
                         baseline_rate=mort_rate)
for (i in 1:length(morthrlist)) {
    dset <- morthrlist[[i]]
    risk <- names(dset)[!names(dset)%in%c('HR', 'stat')]
    names(dset)[!names(dset)%in%c('HR', 'stat')] <- 'Group'
    dset <- data.frame(Risk=c(risk, rep('', nrow(dset)-1)),
                       dset)
    if (i==1) morthrs <- dset else morthrs <- rbind(morthrs, dset)
}
levels(morthrs$Risk) <- list('Clinical Stage'='stage',
                             'Sex'='male')
morthrs$Group <- revalue(as.factor(morthrs$Group),
                         c('0'='Female',
                           '1'='Male'))

 
############################################################
# prep model outputs 
############################################################
# function
makeCI <- function(estimate, lower, upper, fmat) {
    v <- paste0(sprintf(fmat, estimate),
                ' (',
                sprintf(fmat, lower),
                '-',
                sprintf(fmat, upper),
                ')')
    return(v)
}

# survival
surv.in <- read.table('../output/survival_summary.csv',
                      header=TRUE,
                      sep=',')
surv <- data.frame('Measure'=surv.in$Measure,
                   'Mode'=surv.in$Group,
                   'Mean Survival'=makeCI(surv.in$Mean.Survival,
                                 surv.in$Lower,
                                 surv.in$Upper,
                                 '%#.4g'),
                   'Median Survival'=makeCI(surv.in$Median.Survival,
                                   surv.in$Lower.1,
                                   surv.in$Upper.1,
                                   '%#.4g'),
                   'Five Year'=makeCI(surv.in$X5.Year.Survival,
                                      surv.in$Lower.2,
                                      surv.in$Upper.2,
                                      '%#.3g'),
                   'Ten Year'=makeCI(surv.in$X10.Year.Survival,
                                     surv.in$Lower.3,
                                     surv.in$Upper.3,
                                     '%#.3g'),
                   'Twenty Year'=makeCI(surv.in$X20.Year.Survival,
                                        surv.in$Lower.4,
                                        surv.in$Upper.4,
                                        '%#.3g'))
for (i in 3:ncol(surv))
    surv[, i][grepl('NA', surv[, i])] <- NA
surv$Mode <- revalue(surv$Mode,
                     c('clinical'='Clinical',
                       'screened'='Screened'))

# person-years saved
pys.in <- read.table('../output/person_years_saved.csv',
                     header=TRUE,
                     sep=',')
pys.mean <- round(pys.in[1, 2])
pys.lower <- round(pys.in[2, 2])
pys.upper <- round(pys.in[3, 2])
   
# cause-specific mortality
csm <- read.table('../output/net_cause_specific_survival.csv',
                  header=TRUE,
                  sep=',')
csm$Estimate <- makeCI(csm$Rate,
                       csm$X2.5.,
                       csm$X97.5.,
                       '%#.3g')
csm <- subset(csm,
              select=c(Covariate, Value, Mode, Estimate))
names(csm) <- c('Risk', 'Group', 'Mode', 'Estimate')
csm$Risk <- revalue(csm$Risk,
                    c('All'='All Individuals',
                      'stage'='Stage',
                      'male'='Sex'))
csm$Group <- revalue(csm$Group,
                     c('All'='-',
                       '0'='Female',
                       '1'='Male'))
csm$Mode <- revalue(csm$Mode,
                    c('All'='-'))

csm.in <- within(morthrs, 
                 {rate=HR*mort_rate})
csm <- merge(csm, subset(csm.in, select=c(Group, rate)), all=TRUE)
csm <- rename(csm, c('Estimate'='Observed',
                     'rate'='User Input'))
csm <- csm[c(3, 4, 5, 6, 7, 8, 1, 2), ]
csm <- subset(csm, 
              select=c('Risk', 'Group', 'Mode', 'User Input', 'Observed'))
csm$Risk <- as.character(csm$Risk)
csm$Risk[c(2, 4, 5, 6, 8)] <- ''

# delay in cause-specific mortality
delay <- read.table('../output/delay_in_cs_mortality.csv',
                    header=TRUE,
                    sep=',')
delay$Estimate <- makeCI(delay$Mean,
                         delay$X2.5.,
                         delay$X97.5.,
                         '%#.4g')
delay <- subset(delay, select=c(stage, Estimate))
delay <- rename(delay, c(stage='Clinical Stage'))



