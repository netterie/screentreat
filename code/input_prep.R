############################################################
# Prepare inputs for screentreat project
# Breast cancer example
# Current plan: have a wrapper file source this one
# JKB 9/19/2014
############################################################

# Toggle which prep is run
# Options: 'seer_receptor', 'seer_receptorHistStage'
run_this <- 'seer_receptorHistStage'


############################################################
# Receptor frequencies - Historic Stage
############################################################
receptor_file = '~/screentreat/data/bc_2010_receptorfreqHistStage.csv'
receptor_out = '~/screentreat/data/bc_2010_receptorfreq_summaryHistStage.csv'

if (run_this=='seer_receptorHistStage') {

    # Read in file
    dat <- read.csv(receptor_file, header=TRUE, stringsAsFactors=FALSE)
    colnames(dat)

    # Stages to ignore
    ignore <- c('In.situ', 'Localized.regional..Prostate.cases.', 'Unstaged', 'Blank.s.')

    # Original names
    origstage <- c('Localized', 'Regional', 'Distant')

    # Convert 
    # Format
    dat <- dat[,!names(dat)%in%ignore]
    dat <- within(dat, {
                #Early <- S_I
                Early <- Localized
                Advanced <- Regional + Distant
           })
    dat <- dat[,!names(dat)%in%origstage]

    # Finalize for ER+/- only
    er <- subset(dat, !ER_Status%in%c('Unknown', 'Not 1990+ Breast'))
    er <-  within(er, {
            ER_Status[ER_Status%in%c('Negative', 'Borderline')] <- 
                'Negative/Borderline'
        })
    er <- aggregate(formula=cbind(Early, Advanced)~ER_Status, data=er,
                    FUN=sum)
    er$HER2_Status <- NA

    # Finalize for ER+/- and HER2 +/-
    erher2 <- subset(dat, !ER_Status%in%c('Unknown', 'Not 1990+ Breast')
                     & !HER2_Status%in%c('Unknown', 'Not 2010+ Breast'))
    erher2 <-  within(erher2, {
            ER_Status[ER_Status%in%c('Negative', 'Borderline')] <- 
                'Negative/Borderline'
            HER2_Status[HER2_Status%in%c('Negative', 'Borderline')] <- 
                'Negative/Borderline'
        })
    erher2 <- aggregate(formula=cbind(Early, Advanced)~ER_Status+HER2_Status,
                        data=erher2, FUN=sum)

    finaltab <- rbind(er, erher2)[,c('ER_Status', 
                                     'HER2_Status', 'Early', 'Advanced')]
    
    finaltab$Source <- c('SEER 2010. Early=Localized, Advanced=Regional+Distant. Excludes unknown stage, and in each aggregation, unknown receptors.', 
                         'When HER2=NA, numbers reflect marginal over all HER2 including the unknowns',
                         rep('', nrow(finaltab)-2))
    
    write.csv(finaltab,
              file=receptor_out,
              row.names=FALSE)
}
############################################################
# Receptor frequencies
############################################################
receptor_file = '~/screentreat/data/bc_2010_receptorfreq.csv'
receptor_out = '~/screentreat/data/bc_2010_receptorfreq_summary.csv'

if (run_this=='seer_receptor') {

    # Read in file
    dat <- read.csv(receptor_file, header=TRUE, stringsAsFactors=FALSE)
    colnames(dat)

    # Stages to ignore
    ignore <- c('S_0', 'S_NA', 'S_UNK_Stage', 'S_Blank')

    # Original names
    origstage <- grep('S_', names(dat), value=TRUE)

    # Convert 
    # Format
    dat <- dat[,!names(dat)%in%ignore]
    dat <- within(dat, {
                #Early <- S_I
                Early <- S_I + S_IIA
                #Advanced <- S_IIA + S_IIB + S_IIINOS + S_IIIA + S_IIIB + S_IIIC + S_IV
                Advanced <- S_IIB + S_IIINOS + S_IIIA + S_IIIB + S_IIIC + S_IV
           })
    dat <- dat[,!names(dat)%in%origstage]

    # Finalize for ER+/- only
    er <- subset(dat, !ER_Status%in%c('Unknown', 'Not 1990+ Breast'))
    er <-  within(er, {
            ER_Status[ER_Status%in%c('Negative', 'Borderline')] <- 
                'Negative/Borderline'
        })
    er <- aggregate(formula=cbind(Early, Advanced)~ER_Status, data=er,
                    FUN=sum)
    er$HER2_Status <- NA

    # Finalize for ER+/- and HER2 +/-
    erher2 <- subset(dat, !ER_Status%in%c('Unknown', 'Not 1990+ Breast')
                     & !HER2_Status%in%c('Unknown', 'Not 2010+ Breast'))
    erher2 <-  within(erher2, {
            ER_Status[ER_Status%in%c('Negative', 'Borderline')] <- 
                'Negative/Borderline'
            HER2_Status[HER2_Status%in%c('Negative', 'Borderline')] <- 
                'Negative/Borderline'
        })
    erher2 <- aggregate(formula=cbind(Early, Advanced)~ER_Status+HER2_Status,
                        data=erher2, FUN=sum)

    finaltab <- rbind(er, erher2)[,c('ER_Status', 
                                     'HER2_Status', 'Early', 'Advanced')]
    
    finaltab$Source <- c('SEER 2010. Early=AJCC I&IIA, Advanced=AJCC IIB-IV. Excludes unknown stage, and in each aggregation, unknown receptors.', 
                         'When HER2=NA, numbers reflect marginal over all HER2 including the unknowns',
                         rep('', nrow(finaltab)-2))
    
    write.csv(finaltab,
              file=receptor_out,
              row.names=FALSE)
}






# THE FOLLOWING IS OUTDATED (10/27/14):
if (1==0) {
############################################################
# EXAMPLE: receptor_1
############################################################

if (version=='receptor_1') {
    # Groups by which we enter data
    stage <- c('early', 'late')
    subgroup <- c('ER+', 'ER-')
    treatment <- c('None', 'Chemo', 'Tamoxifen', 'Both')

    # Prevalence and baseline rates for stage-subgroup
    stage_subgroup <- expand.grid(subgroup, stage)
    stage_subgroup <- data.frame(stage=stage_subgroup$Var2,
                                 subgroup=stage_subgroup$Var1,
                                 prop=c(0.37, 0.1, 0.39, 0.14),
                                 baseline=c(.021, .021, .102, .102))

    # Prevalence of stage-subgroup-treatment
    stage_subgroup_treat <- expand.grid(treatment, subgroup, stage)
    stage_subgroup_treat <- subset(transform(stage_subgroup_treat,
                                    stage=Var3, 
                                    subgroup=Var2, 
                                    era='old',
                                    treatment=Var1,
                                    prop=NA),
                          select=c('stage', 'subgroup', 'era', 'treatment', 'prop'))
    stage_subgroup_treat <- rbind(stage_subgroup_treat, 
                                  transform(stage_subgroup_treat, era='new'))

    # Hazard ratios for treatment (could always modify these to be for subgroup-treatment or 
    # stage-subgroup-treatment
    treat_HRs <- data.frame(treatment=treatment,
                            HR=c(1, 0.775, 0.71, 0.71*0.775))
                               
}
}
