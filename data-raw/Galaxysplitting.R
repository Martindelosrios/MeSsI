# Reading the data

galaxies      <- read.table('../data/GalaxiesDataset.RData', header = TRUE)
trainClusters <- read.table('../data/ClustersTrainingset.RData', header = TRUE)
testClusters  <- read.table('../data/ClustersTestingset.RData', header = TRUE)
valClusters   <- read.table('../data/ClustersValidationset.RData', header = TRUE)

for(i in 1:length(valClusters$ngroup)){
  Group <- subset(galaxies, galaxies$ngroup == valClusters$ngroup[i])

  if(length(Group$ngroup) == 0){print(valClusters$ngroup[i])}
  if(i == 1){
    valset <- Group
  } else {
    valset <- rbind(valset, Group)
  }
}

for(i in 1:length(testClusters$ngroup)){
  Group <- subset(galaxies, galaxies$ngroup == testClusters$ngroup[i])

  if(length(Group$ngroup) == 0){print(testClusters$ngroup[i])}
  if(i == 1){
    testset <- Group
  } else {
    testset <- rbind(testset, Group)
  }
}

for(i in 1:length(trainClusters$ngroup)){
  Group <- subset(galaxies, galaxies$ngroup == trainClusters$ngroup[i])

  if(length(Group$ngroup) == 0){print(trainClusters$ngroup[i])}
  if(i == 1){
    trainset <- Group
  } else {
    trainset <- rbind(trainset, Group)
  }
}

write.table(trainset, file = '../data/GalaxiesTrainingset.RData', row.names = FALSE, quote = FALSE)
write.table(testset, file = '../data/GalaxiesTestingset.RData', row.names = FALSE, quote = FALSE)
write.table(valset, file = '../data/GalaxiesValidationgset.RData', row.names = FALSE, quote = FALSE)

