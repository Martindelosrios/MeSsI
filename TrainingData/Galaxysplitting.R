# Reading the data

galaxies <- read.table('GalaxiesDataset.dat', header = TRUE)

# Some configuration parameters

ngal  <- length(galaxies$ngroup) # Number of clusters
train <- as.integer(ngal*70/100)
test  <- as.integer(ngal*20/100)
val   <- as.integer(ngal*10/100)

# Let's split the dataset into training, testing and validation sets.

set.seed(288890) # For reproductability reasons

train_index <- sample(1:ngal, size = train)
trainset <- galaxies[train_index,]
testset  <- galaxies[-train_index,]
valset   <- testset[1:val,]
testset  <- testset[-(1:val),]

write.table(trainset, file = 'GalaxiesTrainingset.dat', row.names = FALSE, quote = FALSE)
write.table(testset, file = 'GalaxiesTestingset.dat', row.names = FALSE, quote = FALSE)
write.table(valset, file = 'GalaxiesValidationgset.dat', row.names = FALSE, quote = FALSE)
