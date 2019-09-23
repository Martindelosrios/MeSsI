dat  <- read.table('ClustersTrainingset.dat', header = T) 
test <- read.table('ClustersTestingset.dat', header = T) 


dat        <- dat[,-c(1,2,3,4)]
dat$id_mer <- as.factor(dat$id_mer)

ClustersModel <- train(id_mer~., 
                      data = dat, 
                      method = 'rf', 
                      metric = 'Accuracy') 

dat  <- read.table('GalaxiesTrainingset.dat', header = T) 
test <- read.table('GalaxiesTestingset.dat', header = T) 


dat        <- dat[,-c(1,2,3,4)]
dat$id <- as.factor(dat$id)

GalaxiesModel <- train(id~., 
                      data = dat, 
                      method = 'rf', 
                      metric = 'Accuracy') 

