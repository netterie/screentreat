
###### Inputs #####

sim.size <- 10

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
baseline.hazard.trt1 <- c(.021, .021, .102, .102)

names(baseline.hazard.trt1) <- base.categories

#input a vector giving the RR of each treatment to treatment 1 under each of the base categories


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
base.prop <- c(.3, .3, .2, .2)

empty.prop.matrix <- matrix(nrow=length(base.categories), ncol=length(treatments), dimnames=list(base.categories, treatments))
print(empty.prop.matrix)

# Input the proportion of people in each row getting each treatment.
# For each row of the printed matrix input the proportions getting each treatment
# as elements ofvectors with names: row1, row2, row3, ...
# The sum of each vector needs to be 1
row1 <- c(.2, .8, .1)
row2 <- c(.8, .1, .1)
row3 <- c(.2, 0, .8)
row4 <- c(.2, 0, .8)

row.data <- t(
  data.frame(
    # Input each of the names of the row data as arguments below. Ex: 'row1, row2, row3'
    row1, row2, row3, row4
  ))

prop.matrix <- matrix(data=row.data, nrow=length(base.categories), ncol=length(treatments), dimnames=list(base.categories, treatments))
print(prop.matrix)
# The printed matrix should give the proportion of people in each row getting each treatment



##### Coding functions #####

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


##### Simulation functions #####


generate_simulation_data <- function(n){
  index <- 1:(n)
  sim.data <- data.frame(patient.id=index)
  return(sim.data)
}


##### work in progress #####

sim.cohort <- generate_simulation_data(sim.size)

basestats.sim <- rmultinom(sim.size, 1, base.prop)
basestats.sim <- rmultinom.to.vector(basestats.sim, base.categories)
sim.cohort$base.category <- basestats.sim

sim.cohort$treatment <- sapply(sim.cohort$base.category,
                         FUN=function(x, prop.matrix){
                           p.personal.trt <- prop.matrix[x, ]
                           rmultinom.to.vector(rmultinom(1, 1, p.personal.trt), colnames(prop.matrix))
                         },
                         prop.matrix)



sim.cohort$baseline.hazard <- sapply(sim.cohort$base.category,
                                     FUN=function(x, baseline.hazard.trt1){
                                      baseline.hazard.trt1[x]
                                     },
                                     baseline.hazard.trt1)




