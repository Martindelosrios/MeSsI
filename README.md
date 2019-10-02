# MeSsI Algorithm v2.0

This package contain all the neccesary functions to perform the automatic analysis and classification of merging clusters using the method developed in (https://arxiv.org/abs/1509.02524) by Martín de los Rios, Mariano Domínguez, Dante Paz & Manuel Merchán.

# Prerequisites

This software make an extensive usage of the following packages that must be installed: ```randomForest```, ```nortest```, ```cosmoFns```, ```mclust```, ```e1071```, ```beepr```, ```caret```, ```progress```. 

# Installation

You can install the package directly from your R session using the ```install_github``` function from the ```devtools``` package.

``` R
library('devtools')
install_github('MartindelosRios/MeSsI')
```

# Example

``` R
# Loading the MeSsI library.
library('MeSsI')

# Loading the data
data('GalaxiesDataset')

# Let's see the structure of this dataset
str(GalaxiesDataset)

# As you can see this dataset already have all the properties of the galaxies precomputed.
# We will remove this properties and start with a dataset with only the angular positions (ra, dec), 
#  the redshift (z), the identification of the cluster to which the galaxy belongs (id), the color (color)
#  and the r apparent magnitude (mag).

cat <- GalaxiesDataset[, (1:6)]
colnames(cat)[1] <- 'id'
str(cat)

# Then we just can apply the messi functions to this catalog, optionally given a name to the folder where all the 
#  outputs file will be saved.

messi(cat, folder = 'test')
```

# Authors

Martín de los Rios.
