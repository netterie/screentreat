#############################################################################
## File Name: user_options.r

## File Purpose: 
## Author: Leslie Mallinger
## Date: 8/5/2013 
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
times = '5,10,20'
pop_size = 1000
study_year = 2003

# population characteristics at incidence in absence of screening
#input_method = 'covariate_proportions'    # also 'individual_data'
input_method = 'individual_data'
if (input_method=='individual_data') {
    userdat_file = 'input/input_data.csv'
    create_pop_method = 'simple_bootstrap'  # also 'weighted_bootstrap'
    weighted_bootstrap_table = 'stage,prop\nI,0.2\nII,0.1\nIII,0.3\nIV,0.4'
} else if (input_method=='covariate_proportions') {
    continuous_vars = ''
    categorical_chars1 = 'stage,prop\nI,0.311\nII,0.079\nIII,0.249\nIV,0.361'
    categorical_chars2 = 'agegroup,prop\n55-59,0.428\n60-64,0.306\n65-69,0.178\n70-74,0.088'
    categorical_chars3 = 'male,prop\n0,0.41\n1,0.59'
    categorical_chars4 = ''
    categorical_chars5 = ''
} else stop('Input method must be either individual_data or covariate_proportions')

# population stage distribution in presence of screening
scr_stg_dist = 'order,stage,prop\n1,I,0.500\n2,II,0.070\n3,III,0.213\n4,IV,0.217'

# time to cancer death estimation features
mort_param = 'rate'  # options: 'rate', 'median', 'mean', 'ksurv'
mort_k = NA
mort_value = 0.395
mort_covar1 = 'stage,stat\nI,0.085\nII,0.159\nIII,0.395\nIV,0.795'
mort_covar2 = 'male,HR\n0,0.910\n1,1.078'
mort_covar3 = '' 

# time to other-cause death estimation features
ocd_HR = 1

