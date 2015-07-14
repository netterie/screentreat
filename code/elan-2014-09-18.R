
###### Inputs #####

# file paths to the incidence and life tables
inc.table.file <- '~/screentreat/data/elan_2014-08-29_incidence_table.csv'
life.table.file <- '~/screentreat/data/elan_2014-08-29_life_table.csv'

# simulation size
sim.size <- 1000

# input the screening effect (the RR of being diagnosed with advanced stage cancer)
screen.RR.advanced <- .7

# input a vector containing all treatment options
treatments <- c('NoTreat', 'Tamoxifen', 'Chemo')

# input a vector containing all receptor statuses
receptor <- c('ER+','ER-')

# input a vector containing all stage options (do not change. can only accomodate early and late currently)
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
baseline.hazard.trt1 <- c(.021, .021, .102, .102)

names(baseline.hazard.trt1) <- base.categories



treatment.categories <- vector()
for(i in treatments){
  for(j in receptor){
    treatment.categories <- append(treatment.categories, paste(i,j, sep='.'))
  }
}

empty.hazard.matrix <- matrix(nrow=length(base.categories), ncol=length(treatments), dimnames=list(base.categories, treatments))
print(empty.hazard.matrix)

# Input the proportional hazard of mortality for each col of the printed matrix 
# with respect to the baseline hazard under treatment 1
# as vectors with names: col1, col2, col3, ...
# The first column should be a vector of 1's corresponding to the RR of treatment 1 to treatment 1
col1 <- c(1, 1, 1, 1)
col2 <- c(.83, .83, .83, .83)
col3 <- c(.78, .78, .78, .78)


hazard.matrix <- matrix(nrow=length(base.categories), ncol=length(treatments), dimnames=list(base.categories, treatments),
                        data= # Input a vector containing all the names of the columns of data. ex: 'c(col1, col2, ...)'
                          c(col1, col2, col3)
)

print(base.categories)
# input vector where each element is the proportion 
# of the population going into each of the categories 
# printed on the console
# the sum of the elements of the vector must equal 1
base.prop <- c(.37, .1, .39, .14)

empty.prop.matrix <- matrix(nrow=length(base.categories), ncol=length(treatments), dimnames=list(base.categories, treatments))
print(empty.prop.matrix)

# Input the proportion of people in each row getting each treatment.
# For each row of the printed matrix input the proportions getting each treatment
# as elements ofvectors with names: row1, row2, row3, ...
# The sum of each vector needs to be 1
row1 <- c(.42, .54, .04)
row2 <- c(.37, .05, .58)
row3 <- c(.11, .37, .52)
row4 <- c(.08, .02, .89)
## numbers from Angela Mariotto paper

row.data <- t(
  data.frame(
    # Input each of the names of the row data as arguments below. Ex: 'row1, row2, row3'
    row1, row2, row3, row4
  ))

prop.matrix <- matrix(data=row.data, nrow=length(base.categories), ncol=length(treatments), dimnames=list(base.categories, treatments))
print(prop.matrix)
# The printed matrix should give the proportion of people in each row getting each treatment




#reads in incidence data from file
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


# reads in life table data from file
life.table <- read.csv(life.table.file, stringsAsFactors=FALSE)[ , 1:7]

life.table$Age..years. <- substr(life.table$Age..years., 1, 3)
life.table$Age..years. <- as.numeric(gsub("\\D.*","",life.table$Age..years.))

life.table$Number.surviving.to.age.x <-as.character(life.table$Number.surviving.to.age.x)
life.table$Number.surviving.to.age.x <- gsub(",","",life.table$Number.surviving.to.age.x)
life.table$Number.surviving.to.age.x <- as.numeric(life.table$Number.surviving.to.age.x)

life.table <- subset(life.table, life.table$Age..years. >= 50)

# creates a list of the names of the inputs the simulation will directly call
current.inputs <- list(as.name('sim.size'), 
                       as.name('screen.RR.advanced'),
                       as.name('base.categories'), 
                       as.name('base.prop'), 
                       as.name('prop.matrix'), 
                       as.name('baseline.hazard.trt1'),
                       as.name('hazard.matrix'),
                       as.name('incidence.data'), 
                       as.name('life.table'))



##### save inputs #####

# create a name for the set of inputs (no spaces or characters not accepted as names of objects)
input.name <- 'new.screen.inputs'
assign(input.name, list(sim.size, 
                        screen.RR.advanced,
                        base.categories, 
                        base.prop, prop.matrix, 
                        baseline.hazard.trt1, 
                        hazard.matrix, 
                        incidence.data,
                        life.table))

# this can be used to save the inputs for old.noscreen, new.noscreen, old.screen, and new.screen



##### Coding functions #####


# takes the results from a random multinomial distribution of size=1 and corresponds each result with a category name
rmultinom.to.vector <- function(multi.data, category.names){ 
  for(i in 1:length(category.names)){
    multi.data[i, ] <- multi.data[i, ] * i
  }
  temp <- 1:ncol(multi.data)
  temp <- sapply(temp,
                 FUN=function(x, multi.data){sum(multi.data[, x])},
                 multi.data) 
  temp <- sapply(temp,
                 FUN=function(x, category.names){category.names[x]}, 
                 category.names)
  return(temp)
}



##### internal Simulation functions #####


# creates a dataframe with n patients and gives each patient an id
generate_simulation_data <- function(n){
  index <- 1:(n)
  sim.data <- data.frame(patient.id=index)
  return(sim.data)
}


# generates a clinical diagnosis date for each patient based on incidence data
generate_clinical_diag <- function(data, incidence.data){
  index <- runif(nrow(data), min=0, max(incidence.data$incidence.free.survival))
  clinical.diag <- approx(
    x=incidence.data$incidence.free.survival, 
    y=incidence.data$age,
    xout=index)
  data$clinical.diag <- clinical.diag$y
  data <- transform(data,
                    clinical.diag=ifelse(is.na(clinical.diag), 10000, clinical.diag)) #if not diagnosed by age 87, age of diag is set to 10000
  return(data$clinical.diag)
}


# takes data with a variable for mortality rate named pt.mortality.rate, and gives each patient a survival based on an exponential distribution with their specific mortality rate
generate_survival <- function(data){
  data$survival <- lapply(data$pt.mortality.rate, rexp, n=1)
  return(as.numeric(data$survival))
}


# generates other cause death from life table
generate_ocd <- function(data, life.table){
  index <- runif(nrow(data), min(life.table$Number.surviving.to.age.x), max(life.table$Number.surviving.to.age.x))
  ocd <- approx(
    x=life.table$Number.surviving.to.age.x, 
    y=life.table$Age..years.,
    xout=index)
  data$ocd <- ocd$y
  return(data$ocd)
}



##### simulation function #####


run_simulation <- function(sim.size, 
                           screen.RR.advanced,
                           base.categories, 
                           base.prop, prop.matrix, 
                           baseline.hazard.trt1, 
                           hazard.matrix, 
                           incidence.data,
                           life.table){
  set.seed(149123759)
  
browser()
  late.categories <- vector()
  for(j in receptor){
    late.categories <- append(late.categories, paste('late',j,sep='.'))
  } #creates a list of all base categories with late diagnosis
  
  early.categories <- vector()
  for(j in receptor){
    early.categories <- append(early.categories, paste('early',j,sep='.'))
  } #creates a list of all base categories with early diagnosis
  
  sim.cohort <- generate_simulation_data(sim.size)
  
  # create clinical diagnosis time
  sim.cohort$clinical.diag <- generate_clinical_diag(sim.cohort, incidence.data)
  
  # based on the provided multinomial distribution, this sends each patient into a base category
  basestats.sim <- rmultinom(sim.size, 1, base.prop)
  basestats.sim <- rmultinom.to.vector(basestats.sim, base.categories)
  sim.cohort$base.category <- basestats.sim
  
  # based on screen.RR.advanced, moves some of the late group to the corresponding early group category
  sim.cohort$stage.shift.index <- rbinom(sim.size, 1, 1 - screen.RR.advanced)
  sim.cohort$final.category <- apply(sim.cohort, 1,
                                     FUN=function(x, early.categories, late.categories){
                                       # shift stage but retain receptor status
                                       ifelse(x['stage.shift.index']==1 && x['base.category'] %in% late.categories,
                                              early.categories[which(late.categories == x['base.category'])],
                                              x['base.category']
                                              )
                                       },
                                     early.categories, late.categories)
  
  # for each subject, this chooses a treatment from the multinomial distribution corresponding to their base category
  sim.cohort$treatment <- sapply(sim.cohort$base.category,
                                 FUN=function(x, prop.matrix){
                                   p.personal.trt <- prop.matrix[x, ]
                                   return(rmultinom.to.vector(rmultinom(1, 1, p.personal.trt), colnames(prop.matrix)))
                                 },
                                 prop.matrix)
  
  # gives the baseline hazard under treatment 1 corresponding to each patients base category
  sim.cohort$baseline.hazard <- sapply(sim.cohort$final.category,
                                       FUN=function(x, baseline.hazard.trt1){
                                         baseline.hazard.trt1[x]
                                       },
                                       baseline.hazard.trt1)
  
  # based on each patient's base category and treatment type, gives them a proportional hazard under treatment
  sim.cohort$treatment.effect <- apply(sim.cohort, 1, 
                                       FUN=function(x, hazard.matrix){
                                         return(hazard.matrix[x['final.category'], x['treatment']])
                                       },
                                       hazard.matrix)
  
  # multiplies baseline hazard and proportional hazard to give the patient specific mortality rate
  sim.cohort$pt.mortality.rate <- sim.cohort$baseline.hazard * sim.cohort$treatment.effect
  
  # generates survival from diagnosis for each patient, adds diagnosis age
  sim.cohort$disease.survival <- generate_survival(sim.cohort)
  sim.cohort$survival <- sim.cohort$disease.survival + sim.cohort$clinical.diag
  
  sim.cohort$ocd <- generate_ocd(sim.cohort, life.table)
  
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

# full.results <- data.frame()
# append_results_to_table <- function(full.results, screen.RR.advanced, p.rx.old.early, p.rx.old.advanced, p.rx.new.early, p.rx.new.advanced, rx.effect.early, rx.effect.advanced, old.noscreen, old.screen, new.noscreen, new.screen){
#   temp.results <- data.frame(RR.screen = screen.RR.advanced,
#                              p.old.early = p.rx.old.early,
#                              p.old.late = p.rx.old.advanced,
#                              p.new.early = p.rx.new.early,
#                              p.new.late = p.rx.new.advanced,
#                              rx.eff.early = rx.effect.early,
#                              rx.eff.late = rx.effect.advanced)
#   base.deaths <- sum(old.noscreen$survival < 75 & old.noscreen$event==1)
#   screen.deaths <- sum(old.screen$survival < 75 & old.screen$event==1)
#   treat.deaths <- sum(new.noscreen$survival < 75 & new.noscreen$event==1)
#   treat.screen.deaths <- sum(new.screen$survival < 75 & new.screen$event==1)
#   temp.results$old.MRR <- (base.deaths - screen.deaths)/base.deaths
#   temp.results$new.MRR <- (treat.deaths - treat.screen.deaths)/treat.deaths
#   full.results <- rbind(full.results, temp.results)
#   return(full.results)
# }


##### call simulation #####

stop()
# calls the simulation with current inputs
do.call(run_simulation, current.inputs)


# calls the simulation with a saved list of inputs
#example: do.call(run_simulation, input.to.use )
do.call(run_simulation, new.screen)


old.noscreen <- do.call(run_simulation, old.noscreen.inputs)
old.screen <- do.call(run_simulation, old.screen.inputs)
new.noscreen <- do.call(run_simulation, new.noscreen.inputs)
new.screen <- do.call(run_simulation, new.screen.inputs)



##### Results #####

full.surv.tab <- create_full_surv_table(old.noscreen, old.screen, new.noscreen, new.screen)

