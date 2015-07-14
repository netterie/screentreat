#############################################################################
## File Name: run_file.R

## File Purpose: Main code for CANTRANce area (c) (screening)
## Author: Jeanette Birnbaum
## Date: 10/14/2014

## Additional Comments: 
#############################################################################

############################################################
# Establish model version if this file is not being called 
# by a wrapper
############################################################
# TODO
if (1==0) {
if (!'model_version'%in%ls()) {
    rm(list=ls())
    model_version <- 'ex_erspc'
}
}

############################################################
# Setup
############################################################

source(library_file)
library(cantrance)
set.seed(98103)

############################################################
# Prepare inputs
############################################################

# Identify the trial arms 
trials <- gsub('prop_','',
               grep('prop', colnames(treat_chars), value=TRUE))

# Identify subgroups
subgroups <- as.character(unique(control_notreat$subgroup))

# Prepare numeric and text group IDs:
# First for stage-subgroups
control_notreat <- transform(control_notreat,
                             groupno=1:nrow(control_notreat),
                             groupid=paste(stage, subgroup, 
                                           sep='.'))
treat_chars <- transform(treat_chars,
                         groupid=paste(stage, subgroup, sep='.'))
treat_chars <- merge(treat_chars, 
                     subset(control_notreat, 
                            select=-c(stage, subgroup)),
                     by='groupid',
                     all.x=TRUE, sort=FALSE)

# Next for treatment groups
# Also merge on HRs from treat_chars
treat_ids <- data.frame(groupno=1:length(unique(treat_chars$tx)), 
                        groupid=unique(treat_chars$tx))
warning('Assuming HRs for treatment do not vary across stage-subgroups or trials. If this is not valid, the easiest way to adapt is to specify more treatments, e.g. separate Chemo-Historical from Chemo-Contemporary.')
treat_ids <- merge(treat_ids,
                   subset(treat_chars, groupno==1,
                          select=c(tx, HR)),
                   by.x='groupid',
                   by.y='tx',
                   all.x=TRUE, all.y=FALSE,
                   sort=FALSE)

############################################################
# Prepare population of clinically incident cases
############################################################

# First, simulate the indepenent characteristics given in
# pop_chars. Right now it assumes that each element in the list
# refers to a single variable rather than a joint distribution
# of variables. Very easy to generalize to return the row #
# of the original dataset rather than the value of a single
# variable. It would be easy to instead use Leslie's 
# create_pop_list() function, and that would work with
# a more complex age pattern, too
pop_chars_rows <- lapply(pop_chars, function(x, Npop, Nsim) {
                        sim_multinom(nsims=Npop, 
                                     nreps=Nsim,
                                     probs=x$prop,
                                     names=1:nrow(x))
                        }, pop_size, nsim)

# Age at other-cause death
pop_chars[['age']] <- transform(pop_chars[['age']],
                                    birth_year = study_year-age)

ageOC <- calc_ac_lifespan_pop(popdata=pop_chars,
                                bootrows=pop_chars_rows,
                                results_as_matrix=TRUE, 
                                survHR=ocd_HR) 

# Age at clinical incidence
if (!age_is_ageclin) {

    # Format the incidence data
    inc_data <- format_clinical_incidence(incidence_file)

    # Simulate
    ageclin <- sim_clinical_incidence(popdata=pop_chars,
                               bootrows=pop_chars_rows,
                               incidence=inc_data,
                               results_as_matrix=TRUE)
}

############################################################
# Simulate subgroup and stage
############################################################

# For control arms (for screening arms, stage will change)
control_notreat_rows <- sim_multinom(nsims=pop_size, 
                              nreps=nsim, 
                              probs=control_notreat$prop, 
                              names=1:nrow(control_notreat))

############################################################
# Apply stage shift within subgroups
############################################################

# We can still refer to control_notreat, but now we will 
# have a distribution of rows shifted towards the early 
# stages

# First, identify Early-Advanced stage row pairs within
# subgroups
stage_pairs <- sapply(subgroups, function(x, df) which(df$subgroup==x), 
                      control_notreat)
dimnames(stage_pairs) <- list(control_notreat$stage[stage_pairs[,1]],
                              subgroups)

# Now, for each subgroup, use the HR_advanced stage-shift 
# parameter to determine whether a person gets shifted 
# (1=yes, 0=no). We will generate this for everyone, but we
# will only use it for the advanced stage cases.
shift <- matrix(rbinom(nsim*pop_size, 1, prob=1-HR_advanced),
                nrow=pop_size, ncol=nsim)

# Shift stage
screen_notreat_rows <- control_notreat_rows
for (s in subgroups) {
    adv_cases <- screen_notreat_rows==stage_pairs['Advanced',s]
    screen_notreat_rows[adv_cases & shift==1] <- 
        stage_pairs['Early',s]
    #Check result: table(screen_notreat_rows[adv_cases])
}

############################################################
# Simulate treatment received in control and screening
# arms
############################################################

# For control arms first
# For each trial, simulate treatment (referring to treat_ids)
# by stage-subgroups
control_treatments <- sapply(trials, function(
                          x, 
                          control_notreat_rows,
                          treat_chars,
                          treat_ids,
                          pop_size,
                          nsim) {
                            # Define which trial's proportions to use
                            thisprop <- paste('prop', x, sep='_')
                            treat_results <- 
                                sim_treatment_by_subgroup(treat_ids,
                                                          treat_chars,
                                                          control_notreat_rows,
                                                          thisprop,
                                                          pop_size,
                                                          nsim)
                         }, 
                         control_notreat_rows, 
                         treat_chars,
                         treat_ids,
                         pop_size, 
                         nsim,
                         USE.NAMES=TRUE, simplify=FALSE)

# For screening arms, we need to change treatment only for those
# who got stage shifted, so shift==1 and stage==Advanced. 

# First create indicator of needing to change treatment
shift_treatment <- shift==1 & 
                   control_notreat_rows%in%stage_pairs['Advanced',]

# Now get new treatment assignments for all early-stage cases
screen_treatments <- sapply(trials, function(
                          x, 
                          screen_notreat_rows,
                          treat_chars,
                          treat_ids,
                          pop_size,
                          nsim) {
                            # Define which trial's proportions to use
                            thisprop <- paste('prop', x, sep='_')
                            treat_results <- 
                                sim_treatment_by_subgroup(treat_ids,
                                                          treat_chars,
                                                          screen_notreat_rows,
                                                          thisprop,
                                                          pop_size,
                                                          nsim)
                         }, 
                         screen_notreat_rows, 
                         subset(treat_chars, stage=='Early'),
                         treat_ids,
                         pop_size, 
                         nsim,
                         USE.NAMES=TRUE, simplify=FALSE)

# Replace screen_treatments with control_treatments for non-shifted 
# early-stage cases
for (t in trials) {
    screen_treatments[[t]][!shift_treatment] <- 
        control_treatments[[t]][!shift_treatment]
}


############################################################
# Simulate time to mortality according to stage-subgroup
# and treatment
############################################################

# Baseline mortality rates by stage-subgroup
control_baserate <- return_value_from_id(id=control_notreat_rows, 
                                     df=control_notreat,
                                     value='mortrate')
screen_baserate <- return_value_from_id(id=control_notreat_rows, 
                                     df=control_notreat,
                                     value='mortrate')

# Treatment HRs
control_HRs <- sapply(control_treatments,
                       return_value_from_id, 
                       treat_ids,
                       'HR', 
                       USE.NAMES=TRUE, 
                       simplify=FALSE)
screen_HRs <- sapply(screen_treatments,
                       return_value_from_id, 
                       treat_ids,
                       'HR', 
                       USE.NAMES=TRUE, 
                       simplify=FALSE)

# Final mortality rate
control_rate <- sapply(control_HRs,
                       function(x, rate) { x*rate },
                       control_baserate,
                       USE.NAMES=TRUE, 
                       simplify=FALSE)
screen_rate <- sapply(screen_HRs,
                       function(x, rate) { x*rate },
                       screen_baserate, 
                       USE.NAMES=TRUE,
                       simplify=FALSE)

# Simulate time to cancer death for the 1st control group
ttcd <- matrix(rexp(n=rep(1,pop_size*nsim), rate=control_rate[[1]]),
               nrow=pop_size, ncol=nsim)

# Now for other groups
# For the first control group, we will just get back the same times,
# i.e. ttcd==control_ttcd[[1]]
control_ttcd <- sapply(control_rate,
                       sim_same_qexp,
                       oldtime=ttcd,
                       oldrate=control_rate[[1]],
                       prefix='test',
                       USE.NAMES=TRUE,
                       simplify=FALSE)

