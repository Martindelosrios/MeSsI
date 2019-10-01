#-------------------------------------------------------------------------
#-----------------------------  Example 1  -------------------------------
#-------------------------------------------------------------------------
#
# In this example we will load the preproccesed galaxies and galaxy clusters 
#   training and testing set, train a random forest classifier and apply it
#   to the testing set.
#
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

# Loading the MeSsI package
library('MeSsI')


#--------------------------------------------------------------------------
#---------------------------  Galaxy Clusters  ----------------------------
#--------------------------------------------------------------------------

# Loading the clusters datasets
data('ClustersTrainingset') 
data('ClustersTestingset') 

# Let's remove the id of the group, the angular coordinates and transform
#   the id_mer as a factor.
dat        <- ClustersTrainingset[,-c(1, 2, 3)]
dat$id_mer <- as.factor(dat$id_mer)

# Let's train a random forest
ClustersModel <- train(id_mer~., 
                      data = dat, 
                      method = 'rf', 
                      metric = 'Accuracy') 

# Let's save this model and predict the dynamical status of the galaxy clusters
#   of the testint set.
save(ClustersModel, file = 'ClustersModel_trainingset.Rdata')
pred <- predict(ClustersModel, newdata = ClustersTestingset)

# Let's add the probability of being a merging clusters to the data frame
ClustersTestingset$merProb <- pred$'1' 

# Let's estimate some statistics concerning the performance of the ML model.
ClustersResults <- roc(ClustersTestingset, real.name = 'id_mer', pred.name = 'merProb')

# Let's plot the ROC Curve
plot(ClustersResults$fpr, ClustersResults$tpr, pch = 20, xlab = 'FPR', ylab = 'TPR')

#---------------------------------------------------------------------------
#-----------------------------  Galaxies  ----------------------------------
#---------------------------------------------------------------------------

# Loading the galaxies datasets
data('GalaxiesTrainingset') 
data('GalaxiesTestingset') 


# Let's remove the id of the group, the angular coordinates and transform
#   the id_mer as a factor.
dat       <- GalaxiesTrainingset[,-c(1, 2, 3)]
dat$class <- as.factor(dat$class)


# Let's train a random forest
library(doMC)
registerDoMC(cores = 10)

GalaxiesModel <- train(class~., 
                      data = dat, 
                      method = 'rf', 
                      metric = 'Accuracy') 


# Let's save this model and predict the dynamical status of the galaxy clusters
#   of the testint set.
save(GalaxiesModel, file = 'GalaxiesModel_trainingset.Rdata')
pred <- predict(GalaxiesModel, newdata = GalaxiesTestingset, type = 'prob')

# Let's add the probability of being a merging clusters to the data frame
GalaxiesTestingset$classProb <- pred$'1' 

# Let's estimate some statistics concerning the performance of the ML model.
GalaxiesResults <- roc(GalaxiesTestingset, real.name = 'class', pred.name = 'classProb')

# Let's plot the ROC Curve
plot(GalaxiesResults$fpr, GalaxiesResults$tpr, pch = 20, xlab = 'FPR', ylab = 'TPR')

