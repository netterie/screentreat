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
nsim = 50
times = c(10,13,25)
pop_size = 100000
study_year = 2000

############################################################
# Population features
############################################################



####################for running multiple parameters
#define variables
HR_advanced <- NA
instage_screen_benefit_early <- NA
instage_screen_benefit_advanced <- NA
#save all objects in workspace
space <- ls()
space <- append(space, "space")
for(instage_screen_benefit_advanced in c(.5, .7, .9)){
  for(instage_screen_benefit_early in c(.5, .7, .9)){
    for(HR_advanced in c(.85, .5)){
      rm(list=setdiff(ls(), space)) #removes objects created by previous run
      cat('running settings:Shift: ', HR_advanced,
          ', Early Benefit: ', instage_screen_benefit_early, 
          ', Advanced Benefit: ', instage_screen_benefit_advanced)


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
#HR_advanced = 0.85 #reassigned below

# Within stage treatment benefit
#instage_screen_benefit_early = .9 #reassigned below
#instage_screen_benefit_advanced = .9 #reassigned below

# Add lead time? Default is undefined or FALSE
# If true, add mean lead time in years
lead_time = FALSE
if (lead_time) lt_mean = (40/12)

# Treatment HRs and distributions by subgroup-stage
treat_chars = read.csv(treat_file, header=TRUE, 
                        stringsAsFactors=FALSE)

# Survival distribuion: exponential or weibull?
surv_distr = 'weibull'

# Baseline mortality rates and population proportions by
# subgroup-stages. Subgroup stages specified here must
# match those given in the scrtrt_file
control_notreat = data.frame(stage=c('Early', 'Advanced'),
                             subgroup=c('All', 'All'), 
                             mortshape=c(rep(1.019, 4), rep(0.68, 4)),  ## For Weibull distribution
                             mortscale=c(rep(50.699,4),rep(14.810, 4)), ##
                             prop=c(0.04, 0.38, 0.02, 0.06,
                                    0.06, 0.34, 0.03, 0.07))


############################################################
# Other-cause mortality
############################################################

ocd_HR = 1

############################################################
# Run model
############################################################


      write.table(t(c("Stage_Shift", 
                      "Instage_HR_Early", 
                      "Instage_HR_Advanced")),
                append = TRUE,
                file.path(base_path, model_version, 'output', 
                          'cuminc_mrr_newtable.csv'),
                row.names=FALSE,
                col.names=FALSE,
                sep=",")
      write.table(t(c(HR_advanced,
                    instage_screen_benefit_early,
                    instage_screen_benefit_advanced)),
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