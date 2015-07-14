#############################################################################
## File Name: run_file.r 

## File Purpose: Main code for CANTRANce area (c) (screening)
## Author: Leslie Mallinger
## Date: 5/2/2013
## Edited on: 

## Additional Comments: 
#############################################################################

############################################################
# establish model version if this file is not being called 
# by a wrapper
############################################################
if (!'model_version'%in%ls()) {
    rm(list=ls())
    model_version <- 'ex_erspc'
}


############################################################
# run modes:
# specify whether model is being run locally or for the web 
# and whether R package is up to date
############################################################
app <- TRUE
web <- FALSE
package_updated <- TRUE


############################################################
# set up R
############################################################

# install packages based on run modes above
if (app) {
	if (!"cantrance"%in%rownames(installed.packages())) {
		install.packages("C:/CANTRANce_Projects/R_Shared_Code/cantrance_1.1.zip",
						  repos=NULL,
						  lib=Sys.getenv("R_LIBS_USER"))
	}
	for (p in c("ggplot2", "survival", "cmprsk", 
				"reshape", "plyr", "gridExtra", "msm", "xtable")) {
		if (!p%in%rownames(installed.packages())) {
			install.packages(p, 
							 lib=Sys.getenv("R_LIBS_USER"),
							 repos="http://cran.fhcrc.org/")
		}
	}
}

# load cantrance package and set seed
library(cantrance)
set.seed(98103)

# source and copy files as indicated by run modes above
if (!web & !app) {
	# for local runs using R directly: there is a specific folder structure expected
    if (!package_updated) source('../../Rpackage/branches/main/R/cantrance_library.R')
    setwd(file.path('../examples', model_version))
    file.copy('../../code/run_file.R',
              'input/run_file.R',
              overwrite=TRUE)
} else {
	# local runs using the app, or web runs
    if (!package_updated) source('input/cantrance_library.R')
}


############################################################
# source, format, and confirm validity of  user inputs 
############################################################
cat('Testing validity of user inputs...\n')

if (app) {
    source('./input/diag1022131200.r')
    source('./test_inputs.R')
} else if (web) {
    source('input/params.R')
    source('input/test_inputs.R')
} else {
    source('input/user_options.R')
    source('../../../resources/test_inputs.R') 
}
 

############################################################
# prepare population of clinically incident cases
############################################################
cat('Preparing population of clinically incident cases...\n')

if (input_method=='individual_data') {
    # load population data
    userdat <- read.csv(userdat_file, header=TRUE)
    
    # copy data into model folder
    # this may be turned off if copy_data=FALSE
    if (!"copy_data"%in%ls()) copy_data <- TRUE
    if (copy_data & userdat_file!='input/input_data.csv')
        file.copy(userdat_file,
                  'input/input_data.csv',
                  overwrite=TRUE)

    # ungroup ages, if applicable
    if ('agegroup'%in%names(userdat)) {
        userdat = expand_rangevar(dset=userdat,
                                  rangevar='agegroup',
                                  newvar='age',
                                  splitchar='-')
    }

    # calculate year of birth
    userdat <- transform(userdat, birth_year = study_year - age)

    # bootstrap a population of the desired size from the
    # input data 
    pop <- bootstrap_pop(data=userdat,
                         size=pop_size,
                         how=create_pop_method,
                         covar_table=weighted_bootstrap_table)

    # boostrap rows for nsim datasets
    rows <- bootstrap_rows(data=pop,
                           n_sim=nsim,
                           prefix='sim',
                           return_data=FALSE)

    # extract age at clinical incidence
    age_clin <- bootstrap_var(dset=pop,
                              variable='age',
                              rows=rows)

} else {    ## input_method=='covariate_proportions'
    # convert data frame of continuous covariates to list
    if (sum(is.na(continuous_vars))==0) {
        continuous_vars <- split(continuous_vars,
                                 rownames(continuous_vars))
    } else continuous_vars <- NULL

    # compile continuous and categorical characteristics
    clin_chars <- c(continuous_vars, categorical_chars)

    # expand covariate distributions and calculate year of 
    # birth
    pop <- create_pop_list(char_list=clin_chars,
                           n=pop_size,
                           ref_yr=study_year)

    # bootstrap rows for nsim datasets
    rows <- bootstrap_rows_list(popdata=pop,
                                nrows=pop_size,
                                nsims=nsim,
                                pfx='sim')

    # extract age at clinical incidence
    age_clin <- bootstrap_var(dset=pop[[grep('birth_year', pop)]],
                              variable='age',
                              rows=rows[[grep('birth_year', pop)]])
}


############################################################
# apply stage shift due to screening
############################################################
cat('Applying stage shift due to screening...\n')

# compile inputs on proportion in each stage for clinical
# and screened cases
if (input_method=='individual_data') {
    # compile inputs on proportion in each stage for
    # clinical and screened cases
    stg_dist <- rename(merge(scr_stg_dist,
                             transform(data.frame(table(pop$stage)),
                                       Freq=Freq/sum(Freq)),
                             by.x='stage',
                             by.y='Var1'),
                       c(prop='scr_prop',
                         Freq='clin_prop'))

    # bootstrap rows for screening stage
    # corresponds to 'order' variable in stg_dist
    scr_stage <- bootstrap_shifted_stage(stgs=stg_dist$stage,
                    base_pop=pop$stage,
                    new_cts=stg_dist$scr_prop*pop_size,
                    ordr=stg_dist$order,
                    bootrows=rows,
                    prefix='sims')

    # add screening stage options to population data
    # note that pop here has scr_stage as a 2nd list element,
    # unlike other CANTRANce areas
    pop <- list(pop, data.frame(scr_stage=stg_dist$stage))
    rows <- list(rows, scr_stage)

} else {    ## input_method=='covariate_proportions'
    # compile inputs on proportion in each stage for
    # clinical and screened cases
    stg_dist <- rename(merge(ddply(clin_chars[[grep('stage', clin_chars)]],
                                   .(stage),
                                   summarize,
                                   prop=sum(prop)),
                             scr_stg_dist,
                             by='stage'),
                       c(prop.x='clin_prop',
                         prop.y='scr_prop'))

    # bootstrap rows for screening stage
    scr_stage <- bootstrap_shifted_stage(stgs=stg_dist$stage,
                     base_pop=pop[[grep('stage', pop)]]$stage,
                     new_cts=stg_dist$scr_prop*pop_size,
                     ordr=stg_dist$order,
                     bootrows=rows[[grep('stage', pop)]],
                     prefix='sims')
    
    # add screening stage options to population data
    pop <- c(pop, list(data.frame(scr_stage=stg_dist$stage)))
    rows <- c(rows, list(scr_stage))
}

rx.pop <- data.frame(rx=c('Clinical', 'Screened'))
rx.rows <- replicate(nsim, c(rep(1, pop_size), rep(2, pop_size)))
colnames(rx.rows) <- paste0('sim', 1:nsim)

pop.all <- lapply(1:length(pop),
                      function(x) rbind(pop[[x]], pop[[x]]))
pop.all <- c(pop.all, list(rx.pop))

rows.all <- lapply(1:length(rows),
                       function(x)
                           rbind(rows[[x]], rows[[x]]+nrow(pop[[x]])))
rows.all <- c(rows.all, list(rx.rows))


############################################################
# prepare mortality inputs
############################################################
cat('Preparing inputs for mortality estimation...\n')

# estimate baseline mortality rate based on input value
mort_rate <- get_lambda(param=mort_param,
                        values=mort_value,
                        k=mort_k)

# convert covariate-specific statistics to hazard ratios
mort_HR_list <- convert_HR(lst=mort_covars, 
                           surv_param=mort_param, 
                           surv_k=mort_k,
                           baseline_rate=mort_rate)


############################################################
# estimate time from clinical incidence to cancer death 
# using clinical stage
############################################################
cat('Estimating time to cancer death (clinical stage)...\n')

# apply covariate-specific hazard ratios to get hazard
# ratios for each individual
clin_mort_HR <- calc_HR(data=pop,
                        bootrows=rows,
                        covar_list=mort_HR_list)

# estimate times to cancer death assuming an exponential
# process
clin_ttcd <- sim_piecewiseexp_HRs(baseline_rates=
                                      data.frame(rate=mort_rate,
                                                 times=0),
                                  HRs=clin_mort_HR,
                                  prefix='tfcitcd')


############################################################
# estimate time from clinical incidence to cancer death
# using screened stage
############################################################
cat('Estimating time to cancer death (screened stage)...\n')

# apply covariate-specific hazard ratios to get hazard
# ratios for each individual
scr_mort_HR <- calc_HR(data=pop,
                       bootrows=rows,
                       covar_list=lapply(mort_HR_list,
                           function(x)
                               x <- rename(x,
                                           c(stage='scr_stage'))))

# estimate times to screening cancer death based on quantile
# for clinical case
scr_ttcd <- sim_same_qexp(oldtime=clin_ttcd,
                          oldrate=mort_rate*clin_mort_HR,
                          newrate=mort_rate*scr_mort_HR,
                          prefix='tfsitcd')


############################################################
# estimate time from clinical incidence to other-cause death 
############################################################
cat('Estimating time from clinical incidence to other-cause death...\n')

# calculate age at other-cause death
age_ocd <- calc_ac_lifespan_pop(popdata=pop,
                                bootrows=rows,
                                results_as_matrix=TRUE,
                                survHR=ocd_HR)

# calculate time from clinical incidence to other-cause death
ttod <- calc_ac_lifespan_convert(ocmat=age_ocd,
                                 age=age_clin,
                                 from='age')


############################################################
# summarize mortality results 
############################################################
cat ('Summarizing mortality...\n')

# set times from study start to report survival statistics
x_times <- seq(0, time_max, 0.5)

# determine time from clinical incidence to all-cause death
clin_ttad <- comparevals(mat1=clin_ttcd,
                         mat2=ttod,
                         fun='min',
                         prefix='ttad')
scr_ttad <- comparevals(mat1=scr_ttcd,
                        mat2=ttod,
                        fun='min',
                        prefix='ttad')

# compute all-cause survival
acs <- combinepops(list=list(clin_ttad, scr_ttad),
                   varname='inc',
                   values=c('clinical', 'screened'),
                   df=FALSE)
acs_cdf <- summarize_ecdf(data=acs[[1]], 
                          xvalues=x_times,
                          covar=acs[[2]],
                          covar.order=NULL)

# compute net cause-specific survival
css <- combinepops(list=list(clin_ttcd, scr_ttcd),
                   varname='inc',
                   values=c('clinical', 'screened'),
                   df=FALSE)
net_cdf <- summarize_ecdf(data=css[[1]],
                          xvalues=x_times,
                          covar=css[[2]],
                          covar.order=NULL)

# compute crude cause-specific survival (incorporating
# competing risks)
crude_cdf <- summarize_cmprsk(cancer=css[[1]],
                              other=rbind(ttod, ttod),
                              xvalues=x_times,
                              covar=css[[2]],
                              covar.order=NULL)


############################################################
# plot survival data 
############################################################
cat('Plotting survival data...\n')

# plot and save
ggsave(filename='output/survival_curves.png',
       plot=ggplot_summary(list('All-Cause'=acs_cdf$summ,
                                'Crude'=crude_cdf$summ,
                                'Net'=net_cdf$summ),
                           estimate_column='Mean',
                           covar.order=c('clinical', 'screened')),
       width=10,
       height=6,
       dpi=100)


############################################################
# save survival statistics 
############################################################
cat('Saving survival statistics...\n')

# calculate summary survival statistics
surv_stats <- get_surv_summ(list('All-Cause Survival'=acs_cdf,
                                 'Net Survival'=net_cdf,
                                 'Crude Survival'=crude_cdf),
                            times)

# save csv and html versions
write.csv(surv_stats,
          file='output/survival_summary.csv',
          row.names=FALSE)
print(xtable(surv_stats,
             digits=c(0, 0, 0, 1, 1, 1, 1, 1, 1, rep(3, length(times)*3)),
             caption='Mean, median, and k-year all-cause, crude, and net survival, grouped by mode of detection'),
      type='html',
      file='output/survival_summary.html',
      include.rownames=FALSE,
      caption.placement='top')


############################################################
# calculating person-years saved 
############################################################
cat('Calculating person-years saved...\n')

# calculate person-years saved
pys <- calc_pys(cancer=css[[1]],
                other=rbind(ttod, ttod),
                covar=css[[2]],
                covar.order=NULL,
                obs='clinical',
                calc.time=time_max,
                suppress_negative=FALSE,
                colname='')

# save csv and html versions
write.csv(pys,
          file='output/person_years_saved.csv',
          row.names=TRUE)
print(xtable(pys,
             digits=0,
             caption='Total life-years saved by screening'),
      type='html',
      file='output/person_years_saved.html',
      include.rownames=TRUE,
      include.colnames=FALSE,
      caption.placement='top')


############################################################
# calculate cause-specific survival statistics from simulations
############################################################
cat('Calculating cause-specific survival statistics...\n')

# extract person-time lived post-clinical incidence, cancer 
# death status, and person-time at risk of cancer death
ptime_lived_pi <- ifelse(acs[[1]]<time_max, acs[[1]], time_max)
cd_status <- ifelse(css[[1]]==acs[[1]], 1, 0)
ptime_atrisk_cd <- ifelse(ptime_lived_pi<time_max, ptime_lived_pi, time_max)

# calculate a summary of the appropriate survival statistic
# for incidence-free survival
css_summ <- calc_surv_stat(data=pop.all,
                bootrows=rows.all,
                covar_list=mort_HR_list,
                ptime=ptime_atrisk_cd,
                status=cd_status,
                statistic=mort_param,
                k=mort_k,
                t_units='years',
                stage_shift=TRUE,
                rx_name='Mode')

# save csv and html versions
write.csv(css_summ,
          file='output/net_cause_specific_survival.csv',
          row.names=FALSE)
print(xtable(css_summ,
             digits=c(0, 0, 0, 0, 3, 3, 3, 3),
             caption='FIXME'),
      type='html',
      file='output/net_cause_specific_survival.html',
      include.rownames=FALSE,
      caption.placement='top')


############################################################
# calculate time by which screening delays cancer death 
############################################################
cat('Calculating time by which screening delays cancer death...\n')

# summarize delay in time to cancer death by clinical stage
delay_summ <- calc_summ_by_groups(data=pop,
                                  bootrows=rows,
                                  extra_data=scr_ttcd-clin_ttcd,
                                  extra_data_as_is=TRUE,
                                  extra_data_name='delay',
                                  summ_group='stage',
                                  summ_varname='delay')

# save csv and html versions
write.csv(delay_summ,
          file='output/delay_in_cs_mortality.csv',
          row.names=FALSE)
print(xtable(delay_summ,
             digits=c(0, 0, 2, 2, 2, 2),
             caption='FIXME'),
      type='html',
      file='output/delay_in_cs_mortality.html',
      include.rownames=FALSE,
      caption.placement='top')





