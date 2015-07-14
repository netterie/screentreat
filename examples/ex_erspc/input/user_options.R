#############################################################################
## File Name: user_options.r

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
    continuous_vars = 'varname,mean,sd,min,max\nage,65,7,55,74'
    categorical_chars1 = 'stage,grade,prop\nL,low,0.2\nL,high,0.1\nR,low,0.15\nR,high,0.15\nD,low,0.05\nD,high,0.35'
    categorical_chars2 = 'tx,prop\nRP,.2\nRT,.4\nCM,.4'
    categorical_chars3 = 'male,prop\n0,0\n1,1'
    categorical_chars4 = ''
    categorical_chars5 = ''
} else stop('Input method must be either individual_data or covariate_proportions')

# population stage distribution in presence of screening
scr_stg_dist = 'order,stage,prop\n1,L,0.5\n2,R,0.3\n3,D,0.2'

# time to cancer death estimation features
mort_param = 'ksurv'  # options: 'rate', 'median', 'mean', 'ksurv'
mort_k = 5
mort_value = 0.68
mort_covar1 = 'stage,stat\nL,0.86\nR,0.68\nD,0.32'
mort_covar2 = ''
mort_covar3 = 'tx,HR\nRP,0.65\nRT,0.9\nCM,1'

# time to other-cause death estimation features
ocd_HR = 1

