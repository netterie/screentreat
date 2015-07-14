#############################################################################
## File Name: user_options_screentreat.r

## File Purpose: 
## Author: Leslie Mallinger
## Date: 3/20/2013
## Edited on: 

## Additional Comments: 
#############################################################################

############################################################
# establish model type 
############################################################
mtype = 'screening'

############################################################
# establish user options 
############################################################
# simulation features
nsim = 50
times = '5,10,50'
pop_size = 100
study_year = 2000

# population characteristics at incidence in absence of screening
input_method = 'covariate_proportions'    # also 'individual_data'
if (input_method=='individual_data') {
    userdat_file = 'input/input_data.csv'
    create_pop_method = 'weighted_bootstrap'  # also 'simple_bootstrap'
    weighted_bootstrap_table = 'tx,prop\nCM,0.3\nRP,0.3\nRT,0.4'
} else if (input_method=='covariate_proportions') {
    continuous_vars = 'varname,mean,sd,min,max\nage,50,1,50,50'
    categorical_chars1 = 'stage,tx,tx_new,prop
Early,None,None,0.049
Early,None,Tamoxifen,0.441
Early,None,Chemo,0
Early,Tamoxifen,None,0
Early,Tamoxifen,Tamoxifen,0
Early,Tamoxifen,Chemo,0
Early,Chemo,None,0
Early,Chemo,Tamoxifen,0
Early,Chemo,Chemo,0
Advanced,None,None,0
Advanced,None,Tamoxifen,0
Advanced,None,Chemo,0.408
Advanced,Tamoxifen,None,0
Advanced,Tamoxifen,Tamoxifen,0
Advanced,Tamoxifen,Chemo,0
Advanced,Chemo,None,0
Advanced,Chemo,Tamoxifen,0
Advanced,Chemo,Chemo,0.102'
    categorical_chars2 = 'male,prop\n1,0\n0,1'
    categorical_chars3 = ''
    categorical_chars4 = ''
    categorical_chars5 = ''
} else stop('Input method must be either individual_data or covariate_proportions')

# Is this a comparison of two trials?
two_trials=TRUE

# Is age in the data age at clinical incidence? If not, provide incidence table
age_is_ageclin = FALSE
if (!age_is_ageclin) {
    incidence_file = '~/screentreat/data/bc_clinicical_incidence.csv'
}

# population stage distribution in presence of screening
scr_stg_dist = 'order,stage,prop\n1,Early,0.643\n2,Advanced,0.357'

# time to cancer death estimation features
mort_param = 'rate'  # options: 'rate', 'median', 'mean', 'ksurv'
mort_k = ''
mort_value = 0.021
mort_covar1 = 'stage,stat\nEarly,0.021\nAdvanced,0.102'
mort_covar2 = 'tx,HR\nNone,1\nTamoxifen,1\nChemo,1'
mort_covar3 = 'tx_new,HR\nNone,1\nTamoxifen,1\nChemo,1'

# time to other-cause death estimation features
ocd_HR = 1

