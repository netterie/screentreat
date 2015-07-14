library('reshape')
library('survival') 
library('ggplot2')


##### inputs #############################################################

main <- '~/screentreat'
version <- 'ex_elan'
subversion <- '2014-08-29'

inc.table.file <- '~/screentreat/data/elan_2014-08-29_incidence_table.csv'
life.table.file <- '~/screentreat/data/elan_2014-08-29_life_table.csv'

sim.size              <- 100000

advanced.diag.prop <- .51
screen.RR.advanced <- 0.7

base.mortality.rate.early <- 0.021
base.mortality.rate.advanced <- 0.102

# rx.mortality.RR.early <- .73
# rx.mortality.RR.advanced <- .85

# These can either be single numbers or vectors
rx.effect.early <- c(0.6, 0.7, 0.8, 0.9, 1)
rx.effect.advanced <- c(0.6, 0.7, 0.8, 0.9, 1)
#rx.effect.early <- 1
#rx.effect.advanced <- 1

p.rx.old.early    <- 0
p.rx.new.early    <- .9
p.rx.old.advanced <- .2
p.rx.new.advanced <- 1


##### input data ###################################################

incidence.table <- read.csv(inc.table.file)[1:19, 1:2]
incidence.data <- data.frame(
  age = 0:87, 
  incidence.rate = (approx(incidence.table$age.interpolated, incidence.table$Rate, xout=0:87)$y)/100000)
incidence.data$incidencefree.rate <- 1-incidence.data$incidence.rate
incidence.free.survival <- c()
for(i in 1:88){
  incidence.free.survival[i] <- prod(incidence.data$incidencefree.rate[1:i])
}
incidence.data$incidence.free.survival <- incidence.free.survival
rm(incidence.free.survival)

incidence.data <- subset(incidence.data, incidence.data$age >= 50)



life.table <- read.csv(life.table.file, stringsAsFactors=FALSE)[ , 1:7]

life.table$Age..years. <- as.numeric(gsub("[-][0-9]+","",life.table$Age..years.))
life.table$Age..years.[101] <- 100

life.table$Number.surviving.to.age.x <-as.character(life.table$Number.surviving.to.age.x)
life.table$Number.surviving.to.age.x <- gsub(",","",life.table$Number.surviving.to.age.x)
life.table$Number.surviving.to.age.x <- as.numeric(life.table$Number.surviving.to.age.x)

life.table <- subset(life.table, life.table$Age..years. >= 50)


##### simulation functions #######################################

generate_simulation_data <- function(n){
  index <- 1:(n)
  sim.data <- data.frame(patient.id=index)
  return(sim.data)
}

generate_clinical_diag <- function(data, incidence.data){
  index <- runif(nrow(data), min=0, max(incidence.data$incidence.free.survival))
  clinical.diag <- approx(
    x=incidence.data$incidence.free.survival, 
    y=incidence.data$age,
    xout=index)
  data$clinical.diag <- clinical.diag$y
  data <- transform(data,
                    clinical.diag=ifelse(is.na(clinical.diag), 10000, clinical.diag))
  return(data$clinical.diag)
}

generate_stage_at_diag <- function(data, p.advanced){
  data$stage <- rbinom(nrow(data), 1, p.advanced)
  data <- transform(data, stage = ifelse(stage == 1, 'advanced', 'early'))
  return(data$stage)
}

generate_mortality_rate <- function(data, base.mortality.rate.early, base.mortality.rate.advanced, rx.effect.early, rx.effect.advanced, p.rx.early, p.rx.advanced){
  data$is.rx.early <- rbinom(nrow(data), 1, p.rx.early)
  data$is.rx.advanced <- rbinom(nrow(data), 1, p.rx.advanced)
  data <- transform(data, 
                    is.rx = ifelse(stage=='early', is.rx.early, is.rx.advanced))
  data <- transform(data, 
                    pt.base.mortality.rate = ifelse(stage=='early', base.mortality.rate.early, base.mortality.rate.advanced),
                    pt.rx.effect = ifelse(stage=='early', rx.effect.early, rx.effect.advanced))
  data <- transform(data,
                    pt.mortality.rate = ifelse(is.rx==1, pt.base.mortality.rate*pt.rx.effect, pt.base.mortality.rate))
  return(data$pt.mortality.rate)
}

generate_survival <- function(data){
  data$survival <- lapply(data$pt.mortality, rexp, n=1)
  return(as.numeric(data$survival))
}

generate_ocd <- function(data, life.table){
  index <- runif(nrow(data), min(life.table$Number.surviving.to.age.x), max(life.table$Number.surviving.to.age.x))
  ocd <- approx(
    x=life.table$Number.surviving.to.age.x, 
    y=life.table$Age..years.,
    xout=index)
  data$ocd <- ocd$y
  return(data$ocd)
}


run_simulation <- function(sim.size, 
                           advanced.diag.prop, 
                           base.mortality.rate.early, base.mortality.rate.advanced,
                           rx.effect.early,rx.effect.advanced,
                           p.rx.early,p.rx.advanced,
                           life.table, incidence.data){
#  set.seed(121212)
  set.seed(0126548)
  sim.cohort                <- generate_simulation_data(sim.size)
  sim.cohort$clinical.diag  <- generate_clinical_diag(sim.cohort, incidence.data)
  sim.cohort$stage          <- generate_stage_at_diag(sim.cohort, advanced.diag.prop)
  sim.cohort$pt.mortality   <- generate_mortality_rate(sim.cohort, base.mortality.rate.early, base.mortality.rate.advanced, rx.effect.early, rx.effect.advanced, p.rx.early, p.rx.advanced)
  sim.cohort$disease.surv   <- generate_survival(sim.cohort)
  sim.cohort$survival       <- sim.cohort$disease.surv + sim.cohort$clinical.diag
  sim.cohort$ocd            <- generate_ocd(sim.cohort, life.table)
  sim.cohort                <- transform(sim.cohort,
                                     time=ifelse(survival > ocd, ocd, survival),
                                     event=ifelse(survival > ocd, 0, 1))
  return(sim.cohort)
}


##### tabulation functions #######################################

plot_results <- function(data, title){
  model <- survfit(Surv(data$time, data$event) ~ 1)
  plot(model, xlab='age', ylab='survival', main=title, xlim=c(50,100), ylim=c(.9,1))
}

create_survival_table <- function(data){
  survivaltable <- data.frame(time = seq(from=50, to=75, by=5))
  deaths <- numeric()
  for(i in 1:length(survivaltable$time)){ # calculates number of survival times greater than each t
    deaths[i] <- sum(data$survival < survivaltable$time[i] & data$event==1)
  }
  survivaltable$deaths <- deaths
  return(survivaltable)
}

create_full_surv_table <- function(old.noscreen, old.screen, new.noscreen, new.screen){
  base.survtab <- create_survival_table(old.noscreen)
  screen.survtab <- create_survival_table(old.screen)
  treat.survtab <- create_survival_table(new.noscreen)
  treat.screen.survtab <- create_survival_table(new.screen)
  survtab <- data.frame(
    time = base.survtab$time,
    base = base.survtab$deaths,
    screen = screen.survtab$deaths,
    treat = treat.survtab$deaths,
    treat.screen = treat.screen.survtab$deaths)
  return(survtab)
}

full.results <- data.frame()
append_results_to_table <- function(full.results, screen.RR.advanced, p.rx.old.early, p.rx.old.advanced, p.rx.new.early, p.rx.new.advanced, rx.effect.early, rx.effect.advanced, old.noscreen, old.screen, new.noscreen, new.screen){
  temp.results <- data.frame(RR.screen = screen.RR.advanced,
                             p.old.early = p.rx.old.early,
                             p.old.late = p.rx.old.advanced,
                             p.new.early = p.rx.new.early,
                             p.new.late = p.rx.new.advanced,
                             rx.eff.early = rx.effect.early,
                             rx.eff.late = rx.effect.advanced)
  base.deaths <- sum(old.noscreen$survival < 75 & old.noscreen$event==1)
  screen.deaths <- sum(old.screen$survival < 75 & old.screen$event==1)
  treat.deaths <- sum(new.noscreen$survival < 75 & new.noscreen$event==1)
  treat.screen.deaths <- sum(new.screen$survival < 75 & new.screen$event==1)
  temp.results$old.MRR <- (base.deaths - screen.deaths)/base.deaths
  temp.results$new.MRR <- (treat.deaths - treat.screen.deaths)/treat.deaths
  full.results <- rbind(full.results, temp.results)
  return(full.results)
}

###### run simulation ###########################################

rx.effects <- expand.grid(rx.effect.early, rx.effect.advanced)
colnames(rx.effects) <- c('rx.effect.early', 'rx.effect.advanced')

for (i in 1:nrow(rx.effects)) {
    
    cat ('\n Effect of early and advanced treatments are', 
         rx.effects[i,'rx.effect.early'], 'and', rx.effects[i,'rx.effect.advanced'])

    old.noscreen <- run_simulation(
      sim.size, 
      advanced.diag.prop, 
      base.mortality.rate.early, base.mortality.rate.advanced, 
      rx.effects[i,'rx.effect.early'], rx.effects[i,'rx.effect.advanced'],
      p.rx.old.early, p.rx.old.advanced,
      life.table, incidence.data)

    old.screen <- run_simulation(
      sim.size, 
      advanced.diag.prop*screen.RR.advanced, 
      base.mortality.rate.early, base.mortality.rate.advanced, 
      rx.effects[i,'rx.effect.early'], rx.effects[i,'rx.effect.advanced'],
      p.rx.old.early, p.rx.old.advanced,
      life.table, incidence.data)

    new.noscreen <- run_simulation(
      sim.size, 
      advanced.diag.prop, 
      base.mortality.rate.early, base.mortality.rate.advanced, 
      rx.effects[i,'rx.effect.early'], rx.effects[i,'rx.effect.advanced'],
      p.rx.new.early, p.rx.new.advanced,
      life.table, incidence.data)

    new.screen <- run_simulation(
      sim.size, 
      advanced.diag.prop*screen.RR.advanced, 
      base.mortality.rate.early, base.mortality.rate.advanced, 
      rx.effects[i,'rx.effect.early'], rx.effects[i,'rx.effect.advanced'],
      p.rx.new.early, p.rx.new.advanced,
      life.table, incidence.data)

    survtab <- create_full_surv_table(old.noscreen, old.screen, new.noscreen, new.screen)

    full.results <- append_results_to_table(full.results, screen.RR.advanced, p.rx.old.early, p.rx.old.advanced, p.rx.new.early, p.rx.new.advanced, rx.effects[i,'rx.effect.early'], rx.effects[i,'rx.effect.advanced'], old.noscreen, old.screen, new.noscreen, new.screen)

}


##### output results ################################################################


if (nrow(rx.effects)==1) {

    plot_results(old.noscreen, 'Baseline')
    plot_results(old.screen, 'Screening')
    plot_results(new.noscreen, 'Treatment')
    plot_results(new.screen, 'Treatment and Screening')

} else {

    full.results.plot <- subset(full.results, select=c(rx.eff.early, rx.eff.late, old.MRR, new.MRR))
    full.results.plot <- melt(full.results.plot,
                                id.var=c('rx.eff.early', 'rx.eff.late'),
                                measure.var=c('old.MRR', 'new.MRR'),
                                variable='Trial',
                                value.name='MRR')
    full.results.plot <- transform(full.results.plot,
                                   rx.eff.late=factor(rx.eff.late),
                                   Trial=gsub('.MRR', ' trial', Trial),
                                   MRR=1-value)
    
    plot.MRR <- ggplot(full.results.plot, 
                       aes(x=rx.eff.early,
                           y=MRR,
                           group=rx.eff.late)) +
                geom_point(aes(shape=rx.eff.late)) + 
                facet_grid(.~Trial) + 
                geom_line() +
                ggtitle('MRR at 25-year follow-up')


}

###### save results ######################################################

write.csv(full.results, 
          file.path(main, version, 'output', paste(subversion, 'survival table full.csv', sep='_'))
          , row.names=FALSE)

if (nrow(rx.effects)==1) {

    write.csv(survtab, 
              file.path(main, version, 'output', paste(subversion, 'survival table.csv', sep='_'))
              , row.names=FALSE)

    png(file.path(main, version, 'output', paste(subversion, 'baseline.png', sep='_')))
      plot_results(old.noscreen, 'Baseline')
    dev.off()

    png(file.path(main, version, 'output', paste(subversion, 'treatment effect.png', sep='_')))
      plot_results(new.noscreen, 'Treatment')
    dev.off()

    png(file.path(main, version, 'output', paste(subversion, 'screening effect.png', sep='_')))
      plot_results(old.screen, 'Screening')
    dev.off()

    png(file.path(main, version, 'output', paste(subversion, 'treatment and screening effect.png', sep='_')))
      plot_results(new.screen, 'Treatment and Screening')
    dev.off()

    png(file.path(main, version, 'output', paste(subversion, 'combined results.png', sep='_')))
      par(mfrow=c(2,2))
      plot_results(old.noscreen, 'Baseline')
      plot_results(new.noscreen, 'Treatment')
      plot_results(old.screen, 'Screening')
      plot_results(new.screen, 'Treatment and Screening')
    dev.off()

} else {

    ggsave(file.path(main, version, 'output', paste(subversion, 'MRR.png', sep='_'))
           , plot=plot.MRR)

}




###### work in progress #####


if (1==0) {
    # input a vector containing all treatment options
    treatments <- c('NoTreat', 'Tamoxifen', 'Chemo')
    # input a vector containing all receptor statuses
    receptor <- c('ER+','ER-')
    # input a vector containing all stage options
    stage <- c('early','late')

    base.categories <- vector()
    for(i in stage){
      for(j in receptor){
        base.categories <- append(base.categories, paste(i,j, sep='.'))
      }
    }
    print(base.categories)
    #input a vector with baseline hazards under the first treatment option
    #corresponding to the elements of base.categories (printed on the console)
    baseline.trt1 <- c(.021, .021, .102, .102)


    treatment.categories <- vector()
    for(i in treatments){
      for(j in receptor){
        treatment.categories <- append(treatment.categories, paste(i,j, sep='.'))
      }
    }
}
