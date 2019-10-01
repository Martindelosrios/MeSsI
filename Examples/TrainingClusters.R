library('MeSsI')

data('ClustersTrainingset') 
data('ClustersTestingset') 

dat        <- ClustersTrainingset[,-c(1, 2, 3)]
dat$id_mer <- as.factor(dat$id_mer)

ClustersModel <- train(id_mer~., 
                      data = dat, 
                      method = 'rf', 
                      metric = 'Accuracy') 

data('GalaxiesTrainingset') 
data('GalaxiesTestingset') 


dat    <- dat[,-c(1,2,3,4)]
dat$id <- as.factor(dat$id)


library(doMC)
registerDoMC(cores = 10)

GalaxiesModel <- train(id~., 
                      data = dat, 
                      method = 'rf', 
                      metric = 'Accuracy') 

