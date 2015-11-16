############################################################################
## File Name: user_options.R

## File Purpose: Specify inputs for run_file.R
## Author: Jeanette Birnbaum
## Date: 10/14/2014
## Edited on: 

## Additional Comments: 
#############################################################################

############################################################
# Establish model version if this file is not being called 
# by a wrapper
############################################################
# TODO
if (!'using_wrapper'%in%ls()) {
    warning('Empyting the workspace')
    rm(list=ls())
    model_version <- 'breast_ER-HER2_6Weib'
    rootdir <- '~'
    base_path <- file.path(rootdir, '/screentreat/examples')
}
file.remove(file.path(base_path, 'output', 'cuminc_mrr_newtable.csv'))

############################################################
# Input data files
############################################################

treat_file = file.path(base_path, model_version, 'input', 'input.csv')
incidence_file = file.path(rootdir, '/screentreat/data/bc_1975-1979_incidence.csv')
library_file = file.path(rootdir, '/screentreat/code/screentreat_library.R')
#life_table_file = '~/cantrance/diagnostics/data/life_table_Reed2012.csv'

############################################################
# Simulation features
############################################################
nsim = 100
times = c(10,25)
pop_size = 100000
study_year = 2000

############################################################
# Population features
############################################################



####################for running multiple parameters
#define variables
HR_advanced <- NA
instage_screen_benefit_early <- NA
local_baseline_RR <- NA
weibRRcalc <- function(shape, scale, RR){
  return(scale/(RR^(1/shape)))
}
#instage_screen_benefit_advanced <- NA
#save all objects in workspace
space <- ls()
space <- append(space, "space")
for(instage_screen_benefit_early in c(1)){
  for(HR_advanced in c(.75, .5)){
    for(local_baseline_RR in c(1)){
      rm(list=setdiff(ls(), space)) #removes objects created by previous run
      cat('\nrunning settings:Shift: ', HR_advanced,
          ', Early Benefit: ', instage_screen_benefit_early, 
          ', Base Change Early: ', local_baseline_RR)

pop_chars = 
    list(age=data.frame(age=c(50), prop=c(1)),
         male=data.frame(male=c(0), prop=c(1)))

# Is age in the data age at clinical incidence? 
# If not, provide incidence table
age_is_ageclin = FALSE
if (!age_is_ageclin) {
    inc_table = read.csv(incidence_file, header=TRUE, 
                         stringsAsFactors=FALSE)
}

############################################################
# Screening, treatment and cancer mortality
############################################################

# Stage shift
# HR_advanced = 0.85 #reassigned below

# Within stage treatment benefit
#instage_screen_benefit_early = .9 #reassigned above
instage_screen_benefit_advanced = 1

# Add lead time? Default is undefined or FALSE
# If true, add mean lead time in years
lead_time = FALSE
if (lead_time) lt_mean = (40/12)

# Treatment HRs and distributions by subgroup-stage
treat_chars = read.csv(treat_file, header=TRUE, 
                        stringsAsFactors=FALSE)

# Survival distribuion: exponential or weibull?
surv_distr = 'exponential'

# Baseline mortality rates and population proportions by
# subgroup-stages. Subgroup stages specified here must
# match those given in the scrtrt_file
control_notreat = data.frame(stage = c('Early', 'Advanced'),
                             subgroup = c('All', 'All'), 
                             #mortshape = c(0.948, 0.642),   ## For Weibull distribution
                             #mortscale = c(53.615, 15.69),  ##
                             mortrate=c(.01992, .10693),
                             prop = c(0.496, .504))

#adjust distribution for better baseline survival
# control_notreat[control_notreat$stage=='Early',] <- transform(control_notreat[control_notreat$stage=='Early',],
#                              mortscale = weibRRcalc(mortshape, mortscale, local_baseline_RR)) 

############################################################
# Other-cause mortality
############################################################

ocd_HR = 1

############################################################
# Run model
############################################################


      write.table(t(c("Stage_Shift", 
                      "Instage_HR_Early")),
                append = TRUE,
                file.path(base_path, model_version, 'output', 
                          'cuminc_mrr_newtable.csv'),
                row.names=FALSE,
                col.names=FALSE,
                sep=",")
      write.table(t(c(HR_advanced,
                    instage_screen_benefit_early)),
                  append = TRUE,
                  file.path(base_path, model_version, 'output', 
                            'cuminc_mrr_newtable.csv'),
                  row.names=FALSE,
                  col.names=FALSE,
                  sep=",")
      source(file.path(rootdir, '/screentreat/code/run_file_Sensitivity.R'))
  
    }
  }
}