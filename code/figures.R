############################################################
# Prepare figures for screentreat project
# Breast cancer example
# Current plan: have a wrapper file source this one
# JKB 9/19/2014
############################################################

# Toggle which figure(s) are made
run_this <- 
    #'treatment_freq'
    'hypothetical_results_bar'
    #'hypothetical_results_point'



############################################################
# Treatment frequences
############################################################
treatfreq_file = '~/screentreat/examples/breast_ER-HER2_5/input/input.csv'
treatfreq_out = '~/screentreat/examples/breast_ER-HER2_5/output/treatfreq.pdf'
treatfreq_out2 = '~/screentreat/examples/breast_ER-HER2_5/output/treatfreq_marginal.pdf'

if (run_this=='treatment_freq') {

    library(ggplot2)
    library(reshape2)

    # Read in file
    dat <- read.csv(treatfreq_file, header=TRUE, stringsAsFactors=FALSE)
    colnames(dat)

    # Ignore the +Tras treatments and fix the prop_Perfect
    dat <- subset(dat, !grepl('+Tras', tx) & grepl('HER2-', subgroup))

    # Fix the subgroup
    dat <- within(dat, {
                  subgroup <- gsub('\\+$', '', subgroup)
                  subgroup <- gsub('\\-$', '', subgroup)
                  subgroup <- gsub('HER2', '', subgroup)
            })

    # Melt
    dat <- subset(melt(dat, variable.name='Trial'),
                  Trial!='HR')
    dat <- transform(dat, Trial=gsub('prop_', '', Trial))
    dat <- transform(dat,
                     Trial=factor(Trial, 
                                  levels=c('Historical', 
                                                  'Contemp1999',
                                                  'Perfect'),
                                  labels=c('Historical', 
                                           '1999',
                                           'Perfect')),
                     Treatment=factor(tx, levels=c('None', 'Chemo', 
                                            'Tamoxifen', 
                                            'Tamoxifen+Chemo', 
                                            'AI+Chemo'),
                                   labels=c('None', 'Chemo', 
                                            'Tamoxifen', 
                                            'Tamoxifen+Chemo', 
                                            'AI+Chemo')),
                     Stage=factor(stage, levels=c('Early', 'Advanced'),
                                  labels=c('Early', 'Advanced')),
                     subgroup=factor(subgroup, levels=c('ER+', 'ER-'),
                                     labels=c('ER+', 'ER-')))

    # Black and white 4-color scale
    fcolors <- c(
                 #'#d9d9d9',
                 '#f7f7f7', 
                 '#cccccc', '#969696', '#636363', '#252525')

    # Plot, conditional on stage and ER
    treatfreqpl <- 
        ggplot(data=dat, aes(x=Stage, y=value, fill=Treatment)) + 
        geom_bar(stat="identity") + 
        facet_grid(subgroup~Trial) + theme_bw() + 
        scale_fill_manual(values=fcolors) + 
        scale_y_continuous(name='') + 
        scale_x_discrete(name='') + 
        theme(strip.background=element_rect(colour='white', fill='white'),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())

    ggsave(treatfreq_out, treatfreqpl, width=8, height=5)

    # Plot, unconditional on stage and ER, using proportions from
    # examples/breast_ER-HER2_1/input/user_options.R
    dat <- within(dat, {
                  prop <- value
                  prop[subgroup=='ER+' & stage=='Early'] <- 
                    value[subgroup=='ER+' & stage=='Early']*0.44
                  prop[subgroup=='ER-' & stage=='Early'] <- 
                    value[subgroup=='ER-' & stage=='Early']*0.07
                  prop[subgroup=='ER+' & stage=='Advanced'] <- 
                    value[subgroup=='ER+' & stage=='Advanced']*0.40
                  prop[subgroup=='ER-' & stage=='Advanced'] <- 
                    value[subgroup=='ER-' & stage=='Advanced']*0.11
                     })

    # I don't like this plot so I stopped updating it
    treatfreqpl <- 
        ggplot(data=dat, aes(x=Stage, y=prop, fill=tx)) + 
        geom_bar(stat="identity") + 
        facet_grid(subgroup~Trial) + theme_bw() + 
        scale_fill_manual(values=fcolors)

    ggsave(treatfreq_out2, treatfreqpl, width=8, height=5)

}


############################################################
# Hypothetical results - bar graph
############################################################
#hyp_file = '~/screentreat/examples/breast_hypothetical_3/output/cuminc_mrr_full.csv'
#hyp_file = '~/screentreat/examples/breast_hypothetical_3b/output/cuminc_mrr_full.csv'
hyp_file = '~/screentreat/examples/breast_hypothetical_3half/output/cuminc_mrr_full.csv'
#hyp_out = '~/screentreat/examples/breast_hypothetical_3/output/mrrs_bar.pdf'
#hyp_out = '~/screentreat/examples/breast_hypothetical_3b/output/mrrs.pdf'
hyp_out = '~/screentreat/examples/breast_hypothetical_3half/output/mrrs.pdf'

if (run_this=='hypothetical_results_bar') {

    library(ggplot2)
    library(reshape2)

    # Read in file
    dat <- read.csv(hyp_file, header=TRUE, stringsAsFactors=FALSE)
    colnames(dat)

    # Subset
    dat <- subset(dat, Follow.Up.Year==13 & Measure=='MRR' & 
                        Effect=='Effect of screening, same treatment')

    # Specify HR's for early and late stage, based on trial name
    dat <- transform(dat, 
                     Early_HR=1,
                     Advanced_HR=1,
                     Prop=format(Estimate, digits=2))

    dat <- within(dat, {
                    Early_HR[grepl('e1', Trial)] <- 0.25
                    Early_HR[grepl('e2', Trial)] <- 0.5
                    Early_HR[grepl('e3', Trial)] <- 0.75
                    Advanced_HR[grepl('a1', Trial)] <- 0.25
                    Advanced_HR[grepl('a2', Trial)] <- 0.5
                    Advanced_HR[grepl('a3', Trial)] <- 0.75

        })

    # Black and white 4-color scale
    fcolors <- c('#f7f7f7', '#cccccc', '#969696', '#525252')

    # Plot
    hyppl <- ggplot(dat, ) + 
             geom_bar(aes(x=factor(Early_HR), y=Estimate, 
                          fill=factor(Advanced_HR)),
                      stat='identity', position='dodge') + 
             geom_text(aes(x=factor(Early_HR), y=Estimate, ymax=1.1, 
                           fill=factor(Advanced_HR),
                           label=Prop), 
                       position=position_dodge(width=1),
                       size=2.5,
                       vjust=-1.2, color='grey20') + 
             scale_fill_manual(values=fcolors, 
                               name='HR of advanced stage treatment') +
             scale_y_continuous(name='') + 
             scale_x_discrete(name='HR of early stage treatment') + 
             theme(axis.title.x=element_text(size=10, vjust=-0.75)) + 
             theme_bw() + 
             theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())

    ggsave(hyp_out, hyppl, width=8, height=5)
}



############################################################
# Hypothetical results - point + lines
############################################################
hyp_file = '~/screentreat/examples/breast_hypothetical_3/output/cuminc_mrr_full.csv'
#hyp_file = '~/screentreat/examples/breast_hypothetical_3b/output/cuminc_mrr_full.csv'
#hyp_file = '~/screentreat/examples/breast_hypothetical_3half/output/cuminc_mrr_full.csv'
hyp_out = '~/screentreat/examples/breast_hypothetical_3/output/mrrs.pdf'
#hyp_out = '~/screentreat/examples/breast_hypothetical_3b/output/mrrs.pdf'
#hyp_out = '~/screentreat/examples/breast_hypothetical_3half/output/mrrs.pdf'

if (run_this=='hypothetical_results_point') {

    library(ggplot2)
    library(reshape2)

    # Read in file
    dat <- read.csv(hyp_file, header=TRUE, stringsAsFactors=FALSE)
    colnames(dat)

    # Subset
    dat <- subset(dat, Follow.Up.Year==13 & Measure=='MRR' & 
                        Effect=='Effect of screening, same treatment')

    # Specify HR's for early and late stage, based on trial name
    dat <- transform(dat, 
                     Early_HR=1,
                     Advanced_HR=1,
                     Prop=format(Estimate, digits=2))

    dat <- within(dat, {
                    Early_HR[grepl('e1', Trial)] <- 0.25
                    Early_HR[grepl('e2', Trial)] <- 0.5
                    Early_HR[grepl('e3', Trial)] <- 0.75
                    Advanced_HR[grepl('a1', Trial)] <- 0.25
                    Advanced_HR[grepl('a2', Trial)] <- 0.5
                    Advanced_HR[grepl('a3', Trial)] <- 0.75
        })

    # Black and white 4-color scale
    fcolors <- c('#f7f7f7', '#cccccc', '#969696', '#525252')

    # Plot
    hyppl <- ggplot(dat, ) + 
             geom_point(aes(x=factor(Early_HR), y=Estimate, 
                          group=factor(Advanced_HR),
                          shape=factor(Advanced_HR)),
                        stat='identity', 
                        position=position_dodge(width=0.9)) + 
             geom_errorbar(aes(x=factor(Early_HR), 
                               y=Estimate, 
                               ymin=Estimate, 
                               ymax=Estimate, 
                          group=factor(Advanced_HR)),
                        stat='identity', 
                        position=position_dodge(width=0.9)) + 
             geom_text(aes(x=factor(Early_HR), y=Estimate, ymax=1.1, 
                           fill=factor(Advanced_HR),
                           label=Prop), 
                       position=position_dodge(width=1),
                       size=2.5,
                       vjust=-1.2, color='grey20') + 
             scale_shape_manual(values=c(15:17, 3),
                                name='HR of advanced stage treatment') +
             scale_y_continuous(name='', limits=c(0.80, 1.0)) + 
             scale_x_discrete(name='HR of early stage treatment') + 
             theme(axis.title.x=element_text(size=10, vjust=-0.75)) + 
             theme_bw() + 
             theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())

    ggsave(hyp_out, hyppl, width=8, height=5)
}



