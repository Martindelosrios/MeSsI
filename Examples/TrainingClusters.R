dat  <- read.table('../data/ClustersTrainingset.RData', header = T) 
test <- read.table('../data/ClustersTestingset.RData', header = T) 

dat        <- dat[,-c(1,2,3,4)]
dat$id_mer <- as.factor(dat$id_mer)

ClustersModel <- train(id_mer~., 
                      data = dat, 
                      method = 'rf', 
                      metric = 'Accuracy') 

dat  <- read.table('../data/GalaxiesTrainingset.RData', header = T) 
test <- read.table('../data/GalaxiesTestingset.RData', header = T) 


dat        <- dat[,-c(1,2,3,4)]
dat$id <- as.factor(dat$id)


library(doMC)
registerDoMC(cores = 10)

GalaxiesModel <- train(id~., 
                      data = dat, 
                      method = 'rf', 
                      metric = 'Accuracy') 

