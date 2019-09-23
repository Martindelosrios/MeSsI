# Reading the data

clusters <- read.table('ClustersDataset.dat', header = TRUE)

# Some configuration parameters

nclusters <- length(clusters$ngroup) # Number of clusters
train     <- as.integer(nclusters*70/100)
test      <- as.integer(nclusters*20/100)
val       <- as.integer(nclusters*10/100)

# Let's split the dataset into training, testing and validation sets.

set.seed(288890) # For reproductability reasons

train_index <- sample(1:nclusters, size = train)
trainset <- clusters[train_index,]
testset  <- clusters[-train_index,]
valset   <- testset[1:val,]
testset  <- testset[-(1:val),]

write.table(trainset, file = 'ClustersTrainingset.dat', row.names = FALSE, quote = FALSE)
write.table(testset, file = 'ClustersTestingset.dat', row.names = FALSE, quote = FALSE)
write.table(valset, file = 'ClustersValidationgset.dat', row.names = FALSE, quote = FALSE)
