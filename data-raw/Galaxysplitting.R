# Reading the data

load('FullGalaxiesCatalog.RData')
load('GalaxiesDataset.RData')
load('../data/ClustersTrainingset.RData')
load('../data/ClustersTestingset.RData')
load('../data/ClustersValidationset.RData')

GalaxiesDataset$class <- FullGalaxiesCatalog$class

for(i in 1:length(ClustersValidationset$ngroup)){
  Group <- subset(GalaxiesDataset, GalaxiesDataset$ngroup == ClustersValidationset$ngroup[i])

  if(length(Group$ngroup) == 0){print(ClustersValidationset$ngroup[i])}
  if(i == 1){
    GalaxiesValidationset <- Group
  } else {
    GalaxiesValidationset <- rbind(GalaxiesValidationset, Group)
  }
}

for(i in 1:length(ClustersTestingset$ngroup)){
  Group <- subset(GalaxiesDataset, GalaxiesDataset$ngroup == ClustersTestingset$ngroup[i])

  if(length(Group$ngroup) == 0){print(ClustersTestingset$ngroup[i])}
  if(i == 1){
    GalaxiesTestingset <- Group
  } else {
    GalaxiesTestingset <- rbind(GalaxiesTestingset, Group)
  }
}

for(i in 1:length(ClustersTrainingset$ngroup)){
  Group <- subset(GalaxiesDataset, GalaxiesDataset$ngroup == ClustersTrainingset$ngroup[i])

  if(length(Group$ngroup) == 0){print(ClustersTrainingset$ngroup[i])}
  if(i == 1){
    GalaxiesTrainingset <- Group
  } else {
    GalaxiesTrainingset <- rbind(GalaxiesTrainingset, Group)
  }
}

save(GalaxiesTrainingset, file = '../data/GalaxiesTrainingset.RData')
save(GalaxiesTestingset, file = '../data/GalaxiesTestingset.RData')
save(GalaxiesValidationset, file = '../data/GalaxiesValidationset.RData')

