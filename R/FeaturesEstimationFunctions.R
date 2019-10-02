# DresslerShectmanTest
#{{{
#' DresslerShectmanTest
#' @description This function returns the value of the Dressler-Shectman statistic for a galaxy cluster.
#' @param 
#' cat data frame with a the angular coordinates (ra and dec in radians) and the redshift. The columns must be named 'ra', 'dec' and 'z'
#' @return The Delta value of the Dressler-Shectman test corresponding to the galaxy cluster.
#' @export
#' @examples
#' DresslerShectmanTest(cat)

DresslerShectmanTest <- function(cat){
  n   <- length(cat[,1]) # Number of galaxies in the catalog
  x   <- cat$ra # x coordinate
  y   <- cat$dec # y coordinate
  vel <- cat$z*300000. # velocity

  # Estimation of the peculiar velocity of ngal galaxies
  delta <- 1:n
  for(j in 1:n){
    x0       <- x[j]
    y0       <- y[j]
    r        <- get_distance(x, x0, y, y0, ngal = 10)
    delta[j] <- DresslerShectmanGalaxy(vel,r)
  }

  DELTA <- mean(delta)
  return(DELTA)

}
#}}}
# DresslerShectmanTest2
#{{{
#' DresslerShectmanTest2
#' @description This function returns the value of the Dressler-Shectman statistic for a galaxy cluster with ngal = sqrt(n).
#' @param 
#' cat data frame with a the angular coordinates (ra and dec in radians) and the redshift. The columns must be named 'ra', 'dec' and 'z'.
#' @return The Delta value of the Dressler-Shectman test corresponding to the galaxy cluster.
#' @export
#' @examples
#' DresslerShectmanTest2(cat)

DresslerShectmanTest2 <- function(cat){
  
  x   <- cat$ra # x coordinate 
  y   <- cat$dec # y coordinate
  vel <- cat$z*300000. # velocity
  n   <- length(x) # Number of galaxies in the catalog

  # Estimation of the peculiar velocity of ngal galaxies
  delta <- 1:n
  for(j in 1:n){
    x0       <- x[j]
    y0       <- y[j]
    r        <- get_distance(x, x0, y, y0, ngal = floor(sqrt(n)))
    delta[j] <- DresslerShectmanGalaxy(vel,r)
  }
  DELTA <- mean(delta)
  
  return(DELTA)
}
#}}}
# get_distance
#{{{
#' get_distance
#' @description This function returns the distance from ngal galaxies to a galaxy.
#' @param x Vector with the Right ascensions in radians.
#' @param x0 Right ascension of the central galaxy.
#' @param y Vector with the Declinations in radians.
#' @param y0 Declination of the central galaxy.
#' @param ngal Number of galaxies to be returned
#' @return Vector with the distances to the central galaxy.
#' @export
#' @examples
#' get_distance(x, x0, y, y0, ngal)

get_distance <- function(x, x0, y, y0, ngal){
  n  <- length(x)
  r  <- 1:n
  for(i in 1:n){
    r[i] <- sqrt((((x0-x[i])**2)*cos(y0)*cos(y0))+((y[i]-y0)**2))
  }
  r  <- sort(r,index.return=TRUE)
  ri <- r$ix[2:(ngal+1)]
  return(ri)
}
#}}}
# DresslerShectmanIter
#{{{
#' DresslerShectmanIter
#' @description This function computes the Iterative Dressler-Shectman Test for a galaxy cluster and returns the number of iterations untill convergence.
#' @param 
#' cat data frame with a the angular coordinates (ra and dec in radians) and the redshift. The columns must be named 'ra', 'dec' and 'z'.
#' @return Number of iterations untill convergence.
#' @export
#' @examples
#' DresslerShectmanIter(cat)

DresslerShectmanIter <- function(cat){
  ngal  <- length(cat[,1])
  ngal0 <- ngal/2
  niter <- 30
  corte <- 1000
  c     <- 0

  mat <- cat # We will modify mat at each iteration

  # Estimation of the peculiar velocity of ngal galaxies
  for(i in 1:niter){ 
    if(corte > 0){
      if(ngal > ngal0){
          c     <- c+1
          n     <- length(cat$ra)
          delta <- 1:n
          for(j in 1:n){
            x0       <- cat$ra[j]
            y0       <- cat$dec[j]
            ra       <- cat$ra
            dec      <- cat$dec
            vel      <- cat$z*300000.
            r        <- get_distance(ra, x0, dec, y0, ngal = floor(sqrt(n)))
            delta[j] <- DresslerShectmanGalaxy(vel, r)
          }
          mat   <- data.frame(ra, dec, vel, delta)
          mat   <- subset(mat, mat$delta > 0.5*mean(mat$delta))
          corte <- ngal - length(mat$ra)
          ngal  <- length(mat$ra)
      } else {
          c <- -1
      }  
    }
  }
  return(c)
}
#}}}
# DresslerShectmanGalaxy
#{{{ 
#' DresslerShectmanGalaxy
#' @description This function returns the value of the Dressler-Shectman statistic for an individual galaxy with ngal = sqrt(n).
#' @param vr Vector of velocities of all the galaxies in the cluster.
#' @param r Indixes of the neighbour galaxies to be taking into account in the delta estimation.
#' @return The delta value of the Dressler-Shectman test corresponding to the galaxy.
#' @export
#' @examples
#' DresslerShectmanGalaxy(vr, r)

DresslerShectmanGalaxy <- function(vr, r){
  n       <- length(vr)
  vel_loc <- 1:length(r)
  velmed  <- mean(vr) # Mean velocity of the cluster
  disp    <- sd(vr) # Velocity dispersion of the cluster
 
  velloc   <- mean(vr[r]) # Mean local velocity
  disp_loc <- sd(vr[r])   # Local velocity dispersion
  
  delta <- sqrt(((length(r)+1)/(disp**2))*(((velloc-velmed)**2)+((disp_loc-disp)**2)))
  
  return(delta)
}
#}}}
# DresslerShectmanPval
#{{{
#' DresslerShectmanPval
#' @description This function returns the p-value of the Dressler-Shectman statistic for a galaxy cluster with ngal = 10.
#' @return The p-value of the Dressler-Shectman test corresponding to the galaxy cluster.
#' @param 
#' cat data frame with a the angular coordinates (ra and dec in radians) and the redshift. The columns must be named 'ra', 'dec' and 'z'.
#' @export
#' @examples
#' DresslerShectmanPval(cat)

DresslerShectmanPval <- function(cat){
  n     <- length(cat$z*300000.)
  nran  <- 0
  nmont <- 100
  del   <- DresslerShectmanTest(cat)
  mat   <- cat # we will modify mat in each iteration

  for(i in 1:nmont){
    mat$z   <- sample(mat$z, size = n, replace = FALSE)
    del_mon <- DresslerShectmanTest(mat)
    if(del_mon > del){
      nran <- nran+1
    }
  }
  nran <- nran/nmont
  return(nran)
}
#}}}
# DresslerShectmanPval2
#{{{
#' DresslerShectmanPval2
#' @description This function returns the p-value of the Dressler-Shectman statistic for a galaxy cluster with ngal = floor(sqrt(ngal)).
#' @return The p-value of the Dressler-Shectman test corresponding to the galaxy cluster.
#' @param 
#' cat data frame with a the angular coordinates (ra and dec in radians) and the redshift. The columns must be named 'ra', 'dec' and 'z'.
#' @export
#' @examples
#' DresslerShectmanPval2(cat)

DresslerShectmanPval2 <- function(cat){
  n     <- length(cat$z*300000.)
  nran  <- 0
  nmont <- 100
  del   <- DresslerShectmanTest2(cat)
  mat   <- cat # we will modify mat in each iteration

  for(i in 1:nmont){
    mat$z   <- sample(mat$z, size = n, replace = FALSE)
    del_mon <- DresslerShectmanTest2(mat)
    if(del_mon > del){
      nran <- nran+1
    }
  }
  nran <- nran/nmont
  return(nran)
}
#}}}
# DresslerShectmanIndividual 
#{{{
#' DresslerShectmanIndividual
#' @description This function returns the delta value of the Dressler-Shectman statistic each galaxy of a galaxy cluster with ngal = 10.
#' @return The delta value of the Dressler-Shectman test corresponding to each galaxy of a galaxy cluster.
#' @param 
#' group data frame with a the angular coordinates (ra and dec in radians) and the redshift. The columns must be named 'ra', 'dec' and 'z'.
#' @export
#' @examples
#' DresslerShectmanIndividual(group)

DresslerShectmanIndividual <- function(group){
  n <- length(group$ra)

  # Estimation of the delta statistic for each galaxy 
  delta <- 1:n
  for(j in 1:n){
    r        <- get_distance(group$ra, group$ra[j], group$dec, group$dec[j], ngal  = 10)
    delta[j] <- DresslerShectmanGalaxy(group$z*300000., r)
  }

  return(delta)
}
#}}}
# DresslerShectmanIndividual2 
#{{{
#' DresslerShectmanIndividual2
#' @description This function returns the delta value of the Dressler-Shectman statistic each galaxy of a galaxy cluster with ngal = floor(sqrt(ngal)).
#' @return The delta value of the Dressler-Shectman test corresponding to each galaxy of a galaxy cluster.
#' @param 
#' group data frame with a the angular coordinates (ra and dec in radians) and the redshift. The columns must be named 'ra', 'dec' and 'z'.
#' @export
#' @examples
#' DresslerShectmanIndividual2(group)

DresslerShectmanIndividual2 <- function(group){
  n <- length(group$ra)
  # Estimation of the delta statistic for each galaxy 
  delta <- 1:n
  for(j in 1:n){
    r        <- get_distance(group$ra, group$ra[j], group$dec, group$dec[j], ngal  = floor(sqrt(n)))
    delta[j] <- DresslerShectmanGalaxy(group$z*300000., r)
  }

  return(delta)
}
#}}}
# VelocityDispersionGAP
#{{{
#' VelocityDispersionGAP
#' @description This function estimates the velocity dispersion of a galaxy cluster with the GAP method.
#' @return The GAP velocity dispersion of the galaxy cluster.
#' @param x Vector with the velocities of the galaxies.
#' @export
#' @examples
#' VelocityDispersionGAP(x)

VelocityDispersionGAP <- function(x){
  n <- length(x)-1
  w <- 1:2
  g <- 1:2
  x <- sort(x, decreasing=FALSE)

  for(i in 1:n){
    w[i] <- i*(n+1-i)
    g[i] <- x[i+1]-x[i]
  }

  mul <- g*w
  mul <- sum(mul)
  return(mul)
}

#}}}
# ProjectedDistance
#{{{
#' ProjectedDistance
#' @description This function estimates the projected distances between galaxies needed for the estimation of the virial radius of the galaxy cluster.
#' @param 
#' cat data frame with a the angular coordinates (ra and dec in radians) and the redshift. The columns must be named 'ra', 'dec' and 'z'.
#' @return sum(1/d) where d is a vector with the projected distances.
#' @export
#' @examples
#' ProjectedDistance(cat)

ProjectedDistance <- function(cat){
  id  <- sort(cat$z*300000., decreasing = FALSE, index.return=TRUE)$ix
  vel <- cat$z[id]*300000.
  x   <- cat$ra[id] 
  y   <- cat$dec[id] 
  n   <- length(x)
  r   <- 1:2
  d   <- 0
  
  for(i in 1:n){
    if(i < n){
      for(j in (i+1):n){
        if(x[i] != x[j]){
          if(y[i] != y[j]){
            d    <- d+1
      	    r[d] <- (sqrt((y[i]-y[j])**2+(cos(y[i])*(x[i]-x[j]))**2))
            r[d] <- r[d]*D.A(mean(vel/300000.))
          }
        }
      }
    }
  }
  r <- 1/r
  r <- subset(r,r < 10000000)
  r <- sum(r)
  return(r)
}
#}}}
# ClusterFeatures
#{{{
#' ClusterFeatures
#' @description This function estimates the features for a galaxy cluster.
#' @param 
#' group Data frame with a the angular coordinates (ra and dec in radians), the redshift, the color (g-r), the apparent magnitude (r) and the id of the group for each galaxy. The columns must be named 'ra', 'dec' and 'z', 'color', 'mag' and 'id'.
#' @return Data frame with the features of the galaxy cluster.
#' @export
#' @examples
#' ClusterFeatures(group)

ClusterFeatures <- function(group){
  ngroup       <- group$id[1]
  ra           <- mean(group$ra)
  dec          <- mean(group$dec)
  z            <- mean(group$z)
  ngal         <- length(group$ra)
  color        <- mean(group$color)
  sort_mags    <- sort(group$mag, decreasing = FALSE)
  mag_max      <- sort_mags[1]
  gap_mag      <- sort_mags[2] - sort_mags[1]
  Delta        <- DresslerShectmanTest(group)
  Delta2       <- DresslerShectmanTest2(group)
  pval_ds      <- DresslerShectmanPval(group)
  ind          <- DresslerShectmanIter(group)    
  pval_sw      <- shapiro.test(group$z*300000.)$p.val
  pval_sf      <- sf.test(group$z*300000.)$p.val
  pval_ad      <- ad.test(group$z*300000.)$p.val
  pval_cvm     <- suppressWarnings(cvm.test(group$z*300000.))$p.val
  pval_lillie  <- lillie.test(group$z*300000.)$p.val
  pval_pearson <- pearson.test(group$z*300000.)$p.val

  features <- data.frame(ngroup, ra, dec, z, ngal, color, mag_max, gap_mag, Delta, Delta2, pval_ds, ind, pval_sw, pval_sf, pval_ad, pval_cvm, pval_lillie, pval_pearson)
  return(features)
}
  
#}}}
# get_cluster_features
#{{{
#' get_cluster_features
#' @description This function estimates the features for a catalog of galaxy clusters.
#' @param 
#' dat Data frame with a the angular coordinates (ra and dec in radians), the redshift, the color (g-r), the apparent magnitude (r) and the id of the group for each galaxy of each galaxy cluster. The columns must be named 'ra', 'dec' and 'z', 'color', 'mag' and 'id'.
#' @return Data frame with the features of all the galaxy clusters.
#' @export
#' @examples
#' get_cluster_features(dat)

get_cluster_features <- function(dat, ntotal, name.groups){
  counter  <- 0
  ngal.lim <- 30

  pb <- progress_bar$new(total = floor(ntotal/30))
  for(i in 1:ntotal){
    pb$tick()

    if(length(dat$ra) == 0){break}
    groupid <- dat$id[1] # Id if the group that will be studied
    group   <- subset(dat, dat$id == groupid) # Select the galaxies of the group
    dat     <- subset(dat, dat$id != groupid) # Remove the galaxies of the group from the general catalog
    if(length(group$ra) > ngal.lim){ # Cut in the number of member galaxies
      counter  <- counter+1
      features <- ClusterFeatures(group) # Estimates the features
      if(counter == 1){
        Allfeatures <- features
      } else {
        Allfeatures <- rbind(Allfeatures, features)
      }
    }
  }
  write.table(Allfeatures, file = name.groups, row.names = FALSE)
  return(Allfeatures)
}
#}}}
# ClusterFeatures_new
#{{{
#' ClusterFeatures_new
#' @description This function estimates the features for a galaxy cluster.
#' @param 
#' group Data frame with a the angular coordinates (ra and dec in radians), the redshift, the color (g-r), the apparent magnitude (r) and the id of the group for each galaxy. The columns must be named 'ra', 'dec' and 'z', 'color', 'mag' and 'id'.
#' @return Data frame with the features of the galaxy cluster.
#' @export
#' @examples
#' ClusterFeatures(group)

ClusterFeatures_new <- function(group, featuresFunctions, featuresNames = NULL){

# Clusters general properties
  id  <- group$id[1]
  ra  <- mean(group$ra)
  dec <- mean(group$dec)
  z   <- mean(group$z)

# Features estimation
  nfeat    <- length(featuresFunctions)
  features <- 1:nfeat
  for(i in 1:nfeat){
    features[i] <- featuresFunctions[[i]](group)
  }

  featMAT <- data.frame(id, ra, dec, z)
  featMAT <- cbind(featMAT, as.data.frame(t(features)))
  if(is.null(featuresNames) == FALSE){
    colnames(featMAT) <- c('id', 'ra', 'dec', 'z', featuresNames)
  }
  return(featMAT)
}
  
#}}}
# get_cluster_features_new
#{{{
#' get_cluster_features_new
#' @description This function estimates the features for a catalog of galaxy clusters.
#' @param 
#' dat Data frame with a the angular coordinates (ra and dec in radians), the redshift, the color (g-r), the apparent magnitude (r) and the id of the group for each galaxy of each galaxy cluster. The columns must be named 'ra', 'dec' and 'z', 'color', 'mag' and 'id'.
#' @return Data frame with the features of all the galaxy clusters.
#' @export
#' @examples
#' get_cluster_features(dat)

get_cluster_features_new <- function(dat, ntotal, name.groups){

  featuresFunctions <- list(ngalFunction, colorFunction, mag_maxFunction, gap_maxFunction,
                            DresslerShectmanTest, DresslerShectmanTest2, DresslerShectmanPval,
                            DresslerShectmanIter, shapiro.testGroup, sf.testGroup,
                            ad.testGroup, cvm.testGroup,
                            lillie.testGroup, pearson.testGroup)
  featuresNames <- c('ngal', 'color', 'mag_max', 'gap_max', 'Delta', 'Delta2', 'pval_ds',
                     'ind', 'pval_sw', 'pval_sf', 'pval_ad', 'pval_cvm', 'pval_lillie',
                     'pval_pearson')

  counter  <- 0
  ngal.lim <- 30

  pb <- progress_bar$new(total = floor(ntotal/30))
  for(i in 1:ntotal){
    pb$tick()

    if(length(dat$ra) == 0){break}
    groupid <- dat$id[1] # Id if the group that will be studied
    group   <- subset(dat, dat$id == groupid) # Select the galaxies of the group
    dat     <- subset(dat, dat$id != groupid) # Remove the galaxies of the group from the general catalog
    if(length(group$ra) > ngal.lim){ # Cut in the number of member galaxies
      counter  <- counter+1
      features <- ClusterFeatures_new(group, featuresFunctions, featuresNames) # Estimates the features
      if(counter == 1){
        Allfeatures <- features
      } else {
        Allfeatures <- rbind(Allfeatures, features)
      }
    }
  }
  write.table(Allfeatures, file = name.groups, row.names = FALSE)
  return(Allfeatures)
}
#}}}
# GalaxiesFeatures
#{{{
#' GalaxiesFeatures
#' @description This function estimates the features for the galaxies of a galaxy cluster.
#' @param group_gals Data frame with a the angular coordinates (ra and dec in radians), the redshift, the color (g-r), the apparent magnitude (r) and the id of the group for each galaxy. The columns must be named 'ra', 'dec' and 'z', 'color', 'mag' and 'id'.
#' @param group_data Data frame with a the data corresponding to the galaxy cluster. Must be an object with the same data as the output of ClusterFeatures.
#' @return Data frame with the features of the galaxies.
#' @export
#' @examples
#' GalaxiesFeatures(group_gals, group_data)

GalaxiesFeatures <- function(group_gals, group_data){
  ngals <- length(group_gals$id) # Number of member galaxies of this cluster

  x   <- group_gals$ra
  y   <- group_gals$dec
  vel <- group_gals$z*300000.
  
  delta  <- DresslerShectmanIndividual(group_gals)
  delta2 <- DresslerShectmanIndividual2(group_gals)
  
  pval_sw      <- 1:ngals 
  pval_sf      <- 1:ngals
  pval_ad      <- 1:ngals
  pval_cvm     <- 1:ngals
  pval_lillie  <- 1:ngals
  pval_pearson <- 1:ngals

  for(j in 1:ngals){
    r       <- get_distance(x, x[j], y, y[j], ngal = 10)
    vel_loc <- vel[r]
  
    pval_sw[j]      <- shapiro.test(vel_loc)$p.value
    pval_sf[j]      <- sf.test(vel_loc)$p.value
    pval_ad[j]      <- ad.test(vel_loc)$p.value
    pval_cvm[j]     <- cvm.test(vel_loc)$p.value
    pval_lillie[j]  <- lillie.test(vel_loc)$p.value
    pval_pearson[j] <- pearson.test(vel_loc)$p.value
  }
  
  ngroup      <- group_data$ngroup
  ra          <- x
  dec         <- y
  z           <- vel/300000.
  mag         <- group_gals$mag
  color       <- group_gals$color
  sw_gal      <- pval_sw
  sf_gal      <- pval_sf
  ad_gal      <- pval_ad
  cvm_gal     <- pval_cvm
  lillie_gal  <- pval_lillie
  pearson_gal <- pval_pearson
  Delta       <- replicate(ngals, group_data$Delta)
  pval_ds     <- replicate(ngals, group_data$pval_ds)
  ngal        <- replicate(ngals, group_data$ngal)
  sw_cum      <- replicate(ngals, group_data$pval_sw) 
  sf_cum      <- replicate(ngals, group_data$pval_sf)
  ad_cum      <- replicate(ngals, group_data$pval_ad)
  cvm_cum     <- replicate(ngals, group_data$pval_cvm)
  lillie_cum  <- replicate(ngals, group_data$pval_lillie)
  pearson_cum <- replicate(ngals, group_data$pval_pearson)
  col_cum     <- replicate(ngals, group_data$color) 
  mag_cum     <- replicate(ngals, group_data$mag)
  ind_cum     <- replicate(ngals, group_data$ind)

 features <- data.frame(ngroup, ra, dec, z, mag, color, delta, delta2, sw_gal, sf_gal, ad_gal,
                        cvm_gal, lillie_gal, pearson_gal, Delta, pval_ds, ngal, 
                        sw_cum, sf_cum, ad_cum, cvm_cum, lillie_cum, pearson_cum,
                        col_cum, mag_cum, ind_cum)
  
  return(features)
}
#}}}
# get_galaxies_features
#{{{
#' get_galaxies_features
#' @description This function estimates the features for a catalog of galaxies.
#' @param dat Data frame with a the angular coordinates (ra and dec in radians), the redshift, the color (g-r), the apparent magnitude (r) and the id of the group for each galaxy of each galaxy cluster. The columns must be named 'ra', 'dec' and 'z', 'color', 'mag' and 'id'.
#' @param ClustersData Data frame with a the information of the galaxy cluster. Must have the same data as the output of ClusterFeatures.
#' @return Data frame with the features of all the galaxies of the catalog.
#' @export
#' @examples
#' get_galaxies_features(dat, ClustersData)

get_galaxies_features <- function(dat, ClustersData, name.gal){
  nclusters <- length(ClustersData$ngroup) # Total number of galaxy clusters

  pb <- progress_bar$new(total = nclusters)
  for(i in 1:nclusters){
    pb$tick()

    group_gals    <- subset(dat, dat$id == ClustersData$ngroup[i]) # Id of the cluster that will be studied in this iteration
    group_data    <- ClustersData[i,] # Data of the cluster that will be studied in this iteration
    features <- GalaxiesFeatures(group_gals, group_data) # Estimation of the galaxy features
    if(i == 1){
      AllGalaxyfeatures <- features
    } else {
      AllGalaxyfeatures <- rbind(AllGalaxyfeatures, features)
    }
  }
  write.table(AllGalaxyfeatures, file = name.gal, row.names = FALSE, quote = FALSE)
  return(AllGalaxyfeatures)
}

#}}}
# GalaxiesFeatures_new
#{{{
#' GalaxiesFeatures_new
#' @description This function estimates the features for the galaxies of a galaxy cluster.
#' @param group_gals Data frame with a the angular coordinates (ra and dec in radians), the redshift, the color (g-r), the apparent magnitude (r) and the id of the group for each galaxy. The columns must be named 'ra', 'dec' and 'z', 'color', 'mag' and 'id'.
#' @param group_data Data frame with a the data corresponding to the galaxy cluster. Must be an object with the same data as the output of ClusterFeatures.
#' @return Data frame with the features of the galaxies.
#' @export
#' @examples
#' GalaxiesFeatures(group_gals, group_data)

GalaxiesFeatures_new <- function(group, featuresFunctions, featuresNames){
  ngals <- length(group$id) # Number of member galaxies of this cluster

# General properties of the galaxies
  id    <- group$id
  ra    <- group$ra
  dec   <- group$dec
  z     <- group$z
  mag   <- group$mag
  color <- group$color
  
# Estimating features
  nfeat <- length(featuresFunctions)
  mat   <- matrix(0, ncol = nfeat, nrow = ngals)
  for(i in 1:nfeat){
    mat[,i] <- featuresFunctions[[i]](group) 
  }
  featMAT <- data.frame(id, ra, dec, z, mag, color)
  featMAT <- cbind(featMAT, data.frame(mat))
 
  if(is.null(featuresNames) == FALSE){
    colnames(featMAT)[-(1:6)] <- featuresNames
  }

  return(featMAT)
}
#}}}
# get_galaxies_features_new
#{{{
#' get_galaxies_features_new
#' @description This function estimates the features for a catalog of galaxies.
#' @param dat Data frame with a the angular coordinates (ra and dec in radians), the redshift, the color (g-r), the apparent magnitude (r) and the id of the group for each galaxy of each galaxy cluster. The columns must be named 'ra', 'dec' and 'z', 'color', 'mag' and 'id'.
#' @param ClustersData Data frame with a the information of the galaxy cluster. Must have the same data as the output of ClusterFeatures.
#' @return Data frame with the features of all the galaxies of the catalog.
#' @export
#' @examples
#' get_galaxies_features(dat, ClustersData)

get_galaxies_features_new <- function(dat, ClustersData, name.gal){

  featuresFunctions <- list(DresslerShectmanIndividual, DresslerShectmanIndividual2, 
                         shapiro.testGals, sf.testGals, ad.testGals, cvm.testGals,
                         lillie.testGals, pearson.testGals)
  featuresNames     <- c('delta', 'delta2', 'sw_gal', 'sf_gal', 'ad_gal', 'cvm_gal', 
                         'lillie_gal', 'pearson_gal')
  GfeaturesNames    <- c('Delta', 'pval_ds', 'ngal', 'sw_cum','sf_cum', 'ad_cum', 
                         'cvm_cum', 'lillie_cum', 'pearson_cum', 'col_cum', 
                         'mag_cum', 'ind_cum')

  nclusters <- length(ClustersData$ngroup) # Total number of galaxy clusters

  pb <- progress_bar$new(total = nclusters)
  for(i in 1:nclusters){
    pb$tick()

    group_gals    <- subset(dat, dat$id == ClustersData$ngroup[i]) # Id of the cluster that will be studied in this iteration
    group_data    <- ClustersData[i,] # Data of the cluster that will be studied in this iteration
    features  <- GalaxiesFeatures_new(group_gals, featuresFunctions, featuresNames) # Estimation of the galaxy features
    Gfeatures <- GroupFeatures(group_gals, group_data[,c(9,11,5,13,14,15,16,17,18,6,7,12)], GfeaturesNames)
    features  <- cbind(features, Gfeatures)
    if(i == 1){
      AllGalaxyfeatures <- features
    } else {
      AllGalaxyfeatures <- rbind(AllGalaxyfeatures, features)
    }
  }
  write.table(AllGalaxyfeatures, file = name.gal, row.names = FALSE, quote = FALSE)
  return(AllGalaxyfeatures)
}

#}}}
# get_clusters_classification
#{{{
#' get_clusters_classification
#' @description This function classify galaxies clusters.
#' @param ClustersData Data frame with the features of all the galaxy clusters.
#' @param model Machine Learning model that predicts the dynamical status of a galaxy cluster. Can be the output of Train_cluster_model.
#' @return Vector with the probabilities of being a merging cluster.
#' @export
#' @examples
#' get_clusters_classification(ClustersData, model)

get_clusters_classification <- function(ClustersData, model){

  lim <- 0.3 # Treshol in the probability of being a merging cluster

  model_predictions <- predict(model, newdata = ClustersData, type = 'prob')
  return(model_predictions)
}
#}}}
# get_substructures
#{{{
#' get_substructures
#' @description This function look for the substructures inside merging clusters.
#' @param ClustersData Data frame with the features of all the galaxy clusters.
#' @param GalaxiesData Data frame with the features of all the galaxies of the clusters.
#' @param model Machine Learning model that predicts the probability of a galaxy of being part of a substructure inside a merging cluster. Can be the output of Train_galaxy_model.
#' @return Data frame with the properties of the substructures.
#' @export
#' @examples
#' get_substructures(ClustersData, GalaxiesData, model)

get_substructures <- function(ClustersData, GalaxiesData, model, probLimit, folder){

  classification       <- predict(model, newdata = GalaxiesData, type = 'prob')
  GalaxiesData$relProb <- classification[,1]
  GalaxiesData$subProb <- classification[,2]
  MergingClusters      <- subset(ClustersData$ngroup, ClustersData$merProb > probLimit)

  counter <- 0
  if(length(MergingClusters) > 0){

    pb <- progress_bar$new(total = length(MergingClusters))
    for(i in 1:length(MergingClusters)){
      pb$tick()
  
      group <- subset(GalaxiesData, GalaxiesData$ngroup == MergingClusters[i])
      if(length(group$ra) > 10){
        counter <- counter+1

        ra    <- group$ra
        dec   <- group$dec
        vel   <- group$z*300000.
        mag   <- group$mag
        color <- group$color
        delta <- group$delta
        
        delmin        <- 0
        substructures <- SubstructureIdentification(group, folder)
       
        if(counter == 1){
          AllSubstructures <- substructures
        } else {
          AllSubstructures <- rbind(AllSubstructures, substructures)
        }
      }
    }
  }
  if(exists('AllSubstructures') == FALSE){AllSubstructures$group.id <- -99}
  return(AllSubstructures)
}
#}}}
# messi
#{{{
#' messi
#' @description This function classify the clusters and estimates the merging substructures inside them.
#' @param cat Data frame with the catalog of galaxies. It must have the angular positions in radians (ra, dec), the redshift (z), a flag indicating to which clusters it belong (id), the r apparent magnitude (mag) and a the g-r color (color).
#' @param clusters Boolean indicating if the estimation of the galaxy clusters properties must be done. If FALSE a data frame name ClustersDataset must be load. Defatult TRUE.
#' @param galaxies Boolean indicating if the estimation of the galaxies properties must be done. If FALSE a data frame name GalaxiesDataset must be load. Defatult TRUE.
#' @param classification Boolean indicating if the classification of the galaxy clusters must be done. Default TRUE.
#' @param clustersOutput String indicating the name of the output file that will contain the clusters properties.  Default clustersOutput.dat
#' @param galaxiesOutput String indicating the name of the output file that will contain the galaxies properties. Default galaxiesOutput.dat
#' @param classOutput String indicating the name of the output file that will contain the classification properties. Default ClassOutput.dat
#' @param folder String indicating the name of the folder where the files will be saved. Default 'folder'.
#' @param ClustersML Machine Learning model that will be used for the clusters classification. Default uses the already trained model.
#' @param GalaxiesML Machine Learning model that will be used for the galaxies classification. Default uses the already trained model.
#' @param ntotal Integer indicating the number of galaxy clusters inside the catalog. 0 indicates that this number is not known a priori and must be calculated while estimating the galaxy clusters features. Default 0.
#' @param probLimit Number between 0 and 1, that indicates the treshold that must be used for determining if a cluster is in merger or not. Default 0.3. 
#' @return NULL.
#' @export
#' @examples
#' messi(cat)

messi <- function(cat = -99, 
                  clusters = TRUE,
                  galaxies = TRUE,
                  classification = TRUE, 
                  clustersOutput = 'clustersOutput.dat',
                  galaxiesOutput = 'galaxiesOutput.dat',
                  classOutput    = 'ClassOutput.dat',
                  folder    = 'folder', 
                  ClustersML = 'defaultClustersModel',
                  GalaxiesML = 'defaultGalaxiesModel',
                  #TrainsetClusters = 'trainset_cum.dat',
                  #TrainsetGalaxies = 'trainset_gal.dat',
                  ntotal = 0, 
                  probLimit = 0.3){

  # Configuration and storage folders
#{{{
  if(file.exists(folder) == FALSE){
    system(paste('mkdir', folder, sep = ' '))
  }
  mc <- paste(folder, '/merging_clusters', sep = '') 
  if(file.exists(mc) == FALSE){
    system(paste('mkdir', mc, sep = ' '))
  }
  name.groups  <- paste0(folder, '/', clustersOutput)
  name.gal     <- paste0(folder, '/', galaxiesOutput)
  classOutput  <- paste0(folder, '/', classOutput)
#}}}  
  t1 <- proc.time()

# ---------------------------------------------------------------------------------------
# --------------------------- Galaxy Cluster Features ------------------------------------
# ---------------------------------------------------------------------------------------

  if(clusters == TRUE){
    print('Starting the estimation of the galaxy clusters features')
    if(ntotal == 0){ntotal <- length(cat$ra)}
    ClustersData <- get_cluster_features_new(cat, ntotal, name.groups)
  }

# ---------------------------------------------------------------------------------------
# ------------------------------- Galaxies Features ------------------------------------
# ---------------------------------------------------------------------------------------

  if(galaxies == TRUE){
    print('Starting the estimation of the galaxy features')
    GalaxiesData <- get_galaxies_features_new(cat, ClustersData, name.gal)
  }

# ---------------------------------------------------------------------------------------
# --------------------------- Galaxy Cluster Classification ------------------------------
# ---------------------------------------------------------------------------------------

  if(classification == TRUE){
    print('Starting the classification of the galaxy clusters')
    if(ClustersML == 'defaultClustersModel'){ 
      print('Using Defatult Model')    
      data(ClustersModel)
    }
    Classification       <- get_clusters_classification(ClustersData, ClustersModel)
    ClustersData$relProb <- Classification[,1]
    ClustersData$merProb <- Classification[,2]
    write.table(ClustersData, file = name.groups, row.names = FALSE)

    print('Starting the estimation of the substructures')
    if(GalaxiesML == 'defaultGalaxiesModel'){
      print('Using default Model')
      data(GalaxiesModel)
    }

    Substructures <- get_substructures(ClustersData = ClustersData, 
                                       GalaxiesData = GalaxiesData, 
                                       model = GalaxiesModel,
                                       probLimit = probLimit, folder = folder)

    if(Substructures$group.id[1] != -99){
      names <- c('group.id' ,'FirstSubs', 'SecondSubs', 'rvir1', 
                 'rvir2', 'dvel1', 'dvel2', 'mas1', 'mas2',
                 'par1', 'par2', 'tot', 'racen1', 'racen2', 'deccen1',
                 'deccen2', 'velcen1', 'velcen2')
      write.table(Substructures, file = classOutput, col.names = names, row.names = FALSE, quote = FALSE)
    } else {
      print('There are no merging clusters')
    }
  }

  t2 <- proc.time()
  t  <- (t2[3]-t1[3])
  print(paste('The program delay',toString(t/60),'minutes',sep=' '))
}
#}}}
# SubstructureIdentification
#{{{
#' SubstructureIdentification
#' @description This function performs the identification of the substructures inside a galaxy cluster.
#' @param group Data frame with the galaxies of the cluster. 
#' @return Vector with the properties of the substructures.
#' @export
#' @examples
#' SubstructureIdentification(group)


SubstructureIdentification <- function(group, folder){

  group.id <- paste(folder,'/merging_clusters/',toString(group$ngroup[1]),sep='')
  ngal     <- length(group$ra) 
  mat      <- group

  FirstSubs = SecondSubs = rvir1 = rvir2 = dvel1 = dvel2 = mas1 = mas2 <- -99 # Initializate
  par1 = par2 = tot = racen1 = racen2 = deccen1 = deccen2 = velcen1 = velcen2 <- -99
  if(ngal > 10){
    GroupSubs <- WeightedMclust(group) # Mclust pesado por peso
    nSubs     <- max(GroupSubs$SubId)

    if(nSubs > 1){
      name <- paste(toString(group.id),'_galaxies.dat',sep='')
      if(file.exists(name) == FALSE){
        write.table(GroupSubs, file = name, row.names = FALSE, quote = FALSE)
      }

      SubsNgal <- 1:nSubs
      for(j in 1:nSubs){
        SubsNgal[j] <- length(which(GroupSubs$SubId == j))
      }
      SubsNgal   <- sort(SubsNgal, decreasing = TRUE, index.return = TRUE)
      FirstSubs  <- SubsNgal$x[1]
      SecondSubs <- SubsNgal$x[2]
      tot        <- (FirstSubs+SecondSubs)/ngal

      if(SecondSubs > 4){
        group1   <- subset(GroupSubs, GroupSubs$SubId == SubsNgal$ix[1])
        gap_dvel <- VelocityDispersionGAP(group1$z*300000.)
        ProjDist <- ProjectedDistance(data.frame(group1$ra, group1$dec, group1$z*300000.))
        rvir1    <- (pi*(length(group1$ra)-1)*length(group1$ra))/(2*ProjDist)
        dvel1    <- (sqrt(pi)*gap_dvel)/(length(group1$ra)*(length(group1$ra)-1))
        mas1     <- (5*(dvel1**2)*rvir1)/(4.314465e-09)
        racen1   <- mean(group1$ra)*pi/180
        deccen1  <- mean(group1$dec)*pi/180
        velcen1  <- mean(group1$z*300000.)

        group2   <- subset(GroupSubs, GroupSubs$SubId == SubsNgal$ix[2])
        gap_dvel <- VelocityDispersionGAP(group2$z*300000.)
        ProjDist <- ProjectedDistance(data.frame(group2$ra, group2$dec, group2$z*300000.))
        rvir2    <- (pi*(length(group2$ra)-1)*length(group2$ra))/(2*ProjDist)
        dvel2    <- (sqrt(pi)*gap_dvel)/(length(group2$ra)*(length(group2$ra)-1))
        mas2     <- (5*(dvel2**2)*rvir2)/(4.314465e-09)
        racen2   <- mean(group2$ra)*pi/180
        deccen2  <- mean(group2$dec)*pi/180
        velcen2  <- mean(group2$z*300000.)

        deccen    <- (deccen1+deccen2)/2
        velcen    <- (velcen1+velcen2)/2
        dist12    <- (sqrt((deccen1-deccen2)**2+(cos(deccen)*(racen1-racen2))**2))*(velcen/100)/(1*(1+(velcen/300000.)))
        par1      <- dist12/(rvir1+rvir2)
        par2      <- abs(velcen1-velcen2)/(dvel1+dvel2)
      }
    }
  }

  SubsProperties <- c(group.id = group$ngroup[1], FirstSubs, SecondSubs, rvir1, rvir2, dvel1, dvel2, mas1, mas2,
         par1, par2, tot, racen1, racen2, deccen1, deccen2, velcen1, velcen2)
  return(as.data.frame(t(SubsProperties)))
}
#}}}
# WeightedMclust
#{{{
#' WeightedMclust
#' @description This function performs a wegihted mixture of gaussians.
#' @param group Data frame with the galaxies of the cluster. 
#' @return data frame with the same info that the input data frame plus the substructure id.
#' @export
#' @examples
#' WeightedMclust(cat)


WeightedMclust <- function(group){
  for(k in 1:length(group$ra)){
    dmin    <- min(group$subProb)
    dmax    <- max(group$subProb)
    npoints <- floor(((group$subProb[k]-dmin)/(dmax-dmin))*49+1)
    if(k == 1){
      xnew  <- replicate(n = npoints, group$ra[k])
      ynew  <- replicate(n = npoints, group$dec[k])
      GlxId <- sample(k, size = npoints, replace = TRUE)
    } else {
      xnew  <- c(xnew, replicate(n = npoints, group$ra[k]))
      ynew  <- c(ynew, replicate(n = npoints, group$dec[k]))
      GlxId <- c(GlxId, replicate(n = npoints, k))
    }
  }

  auxGroup     <- data.frame(xnew, ynew) 
  mclust_model <- Mclust(auxGroup, G = 1:2)

  auxGroup <- data.frame(SubId = mclust_model$classification, GlxId)
  auxGroup <- unique(auxGroup) # Remove all the replications
  group    <- cbind(group, SubId = auxGroup$SubId)
  return(group)
}
#}}}
# roc
#{{{
#' roc
#' @description This function estimates some statistics for measuring the performance of the machine learning model.
#' @param dat Data frame with the real class and the prediction of the probability.
#' @param real.name String with the name of the variable that is the real class. Defatult is 'class'. 
#' @param pred.name String with the name of the variable that is the predicted probability. Defaults is 'merProb'.
#' @return data frame with the statistics. 
#' @export
#' @examples
#' roc(cat)

roc <- function(dat, real.name = 'class', pred.name = 'merProb'){
  real_merger <- length(which(dat[real.name] == 1))
  real_relax  <- length(which(dat[real.name] == 0))
  treshold    <- seq(from = 0.0001, to = 1, length.out = 10)

  fpr <- 1:length(treshold)
  tpr <- 1:length(treshold)
  tnr <- 1:length(treshold)
  ppv <- 1:length(treshold)
  acc <- 1:length(treshold)
  f1  <- 1:length(treshold)
  for(i in 1:length(treshold)){
    mer <- subset(dat, dat[pred.name] >= treshold[i])
    rel <- subset(dat, dat[pred.name] < treshold[i])

    if(length(mer$ra) > 0){
      pred_mer  <- length(mer[pred.name])
      true_mer  <- length(which(mer[real.name] == 1)) 
      false_mer <- length(which(mer[real.name] == 0)) 
      ppv[i]    <- true_mer/pred_mer
      f1[i]     <- 2*ppv[i]*tpr[i]/(ppv[i]+tpr[i])
    } else {
      pred_mer  <- 0 
      true_mer  <- 0 
      false_mer <- 0 
      ppv[i]    <- -99
      f1[i]     <- -99
    }
    if(length(rel$ra) > 0){
      pred_rel  <- length(rel[pred.name])
      true_rel  <- length(which(rel[real.name] == 0)) 
      false_rel <- length(which(rel[real.name] == 1)) 
    } else {
      pred_rel  <- 0 
      true_rel  <- 0 
      false_rel <- 0 
    } 
    fpr[i] <- false_mer/real_relax
    tpr[i] <- true_mer/real_merger
    tnr[i] <- true_rel/real_relax
    acc[i] <- (true_mer+true_rel)/(real_merger+real_relax)
  }
 
  return(data.frame(treshold, fpr, tpr, tnr, ppv, acc, f1))
}


#}}}
# ngalFunction
#{{{
ngalFunction <- function(group){
  return(length(group$ra))  
}
#}}}
# colorFunction
#{{{
colorFunction <- function(group){
  return(mean(group$color))  
}
#}}}
# mag_maxFunction
#{{{
mag_maxFunction <- function(group){
  return(max(group$mag))
}
#}}}
# gap_maxFunction
#{{{
gap_maxFunction <- function(group){
  sort_mag <- sort(group$mag, decreasing = TRUE, index.return = FALSE)
  return((sort_mag[1]-sort_mag[2]))
}
#}}}
# shapiro.testGroup
#{{{
shapiro.testGroup <- function(group){
  return(shapiro.test(group$z*300000.)$p.val)
}
#}}}
# sf.testGroup
#{{{
sf.testGroup <- function(group){
  return(sf.test(group$z*300000.)$p.val)
}
#}}}
# ad.testGroup
#{{{
ad.testGroup <- function(group){
  return(ad.test(group$z*300000.)$p.val)
}
#}}}
# cvm.testGroup
#{{{
cvm.testGroup <- function(group){
  return(cvm.test(group$z*300000.)$p.val)
}
#}}}
# lillie.testGroup
#{{{
lillie.testGroup <- function(group){
  return(lillie.test(group$z*300000.)$p.val)
}
#}}}
# pearson.testGroup
#{{{
pearson.testGroup <- function(group){
  return(pearson.test(group$z*300000.)$p.val)
}
#}}}
# shapiro.testGals
#{{{
shapiro.testGals <- function(group){

  ngals   <- length(group$ra)
  pval_sw <- 1:ngals 
  for(j in 1:ngals){
    r          <- get_distance(group$ra, group$ra[j], group$dec, group$dec[j], ngal = 10)
    vel_loc    <- group$z[r]*300000.
    pval_sw[j] <- shapiro.test(vel_loc)$p.value
  }
  return(pval_sw)
}
#}}}
# sf.testGals
#{{{
sf.testGals <- function(group){

  ngals   <- length(group$ra)
  pval_sf <- 1:ngals 
  for(j in 1:ngals){
    r       <- get_distance(group$ra, group$ra[j], group$dec, group$dec[j], ngal = 10)
    vel_loc <- group$z[r]*300000.
  
    pval_sf[j]      <- sf.test(vel_loc)$p.value
  }
  return(pval_sf)
}
#}}}
# ad.testGals
#{{{
ad.testGals <- function(group){

  ngals   <- length(group$ra)
  pval_ad <- 1:ngals 
  for(j in 1:ngals){
    r       <- get_distance(group$ra, group$ra[j], group$dec, group$dec[j], ngal = 10)
    vel_loc <- group$z[r]*300000.
  
    pval_ad[j]      <- ad.test(vel_loc)$p.value
  }
  return(pval_ad)
}
#}}}
# cvm.testGals
#{{{
cvm.testGals <- function(group){

  ngals   <- length(group$ra)
  pval_cvm <- 1:ngals 
  for(j in 1:ngals){
    r       <- get_distance(group$ra, group$ra[j], group$dec, group$dec[j], ngal = 10)
    vel_loc <- group$z[r]*300000.
  
    pval_cvm[j]     <- cvm.test(vel_loc)$p.value
  }
  return(pval_cvm)
}
#}}}
# lillie.testGals
#{{{
lillie.testGals <- function(group){

  ngals   <- length(group$ra)
  pval_lillie <- 1:ngals 
  for(j in 1:ngals){
    r       <- get_distance(group$ra, group$ra[j], group$dec, group$dec[j], ngal = 10)
    vel_loc <- group$z[r]*300000.
  
    pval_lillie[j]  <- lillie.test(vel_loc)$p.value
  }
  return(pval_lillie)
}
#}}}
# pearson.testGals
#{{{
pearson.testGals <- function(group){

  ngals   <- length(group$ra)
  pval_pearson <- 1:ngals 
  for(j in 1:ngals){
    r       <- get_distance(group$ra, group$ra[j], group$dec, group$dec[j], ngal = 10)
    vel_loc <- group$z[r]*300000.
  
    pval_pearson[j] <- pearson.test(vel_loc)$p.value
  }
  return(pval_pearson)
}
#}}}
# GroupFeatures 
#{{{
GroupFeatures <- function(group_gals, group_data, featuresNames){
  ngals <- length(group_gals$ra)
  nfeat <- length(group_data)

  mat <- matrix(0, ncol = nfeat, nrow = ngals)
  for(i in 1:nfeat){
    mat[,i] <- replicate(ngals, group_data[,i])
  } 
  mat <- data.frame(mat)
  if(is.null(featuresNames) == FALSE){
    colnames(mat) <- featuresNames
  }
  return(mat)
}
#}}}



#----------------------------------------------------------------
#-------------------------  TO DO  --------------------------
#----------------------------------------------------------------

# Train_cluster_model
#{{{
#}}}
# Train_galaxy_model
#{{{
#}}}

    # RANKING-RELAXED
  #{{{
##  trainset<-read.table(file=name_trainset_cum,header=TRUE)
##  nrank=nrank
##  
##  dat_cum<-read.table(file=name.groups,header=TRUE)
##  if(length(dat_cum$rf.pred.rel) > 0){ #Si ya existe la medicion aca la elimino para luego reemplazarla
##    dat_cum <- subset(dat_cum, select = -c(rf.pred.rel))
##  }
##  lim=0.6
##  
##  # Lectura de datos de cumulos
##  #if(mean(dat_cum$z) > 0.15){
##  #  trainset_rel<-read.table('/media/martin/store1/trabajos/mock/guo/mock_cumulos_snap56.dat',header=TRUE)
##  #} else {
##  #  trainset_rel<-read.table('/media/martin/store1/trabajos/mock/guo/mock_cumulos_snap63.dat',header=TRUE)
##  #}
##  cont.rel=0
##  rf.pred.aux<-1:length(dat_cum$ngroup)
##  rf.pred.aux[]=0
##  
##  pb <- txtProgressBar(title = "progress bar", min = 0,max = nrank, width = 82)
##  for(jj in 1:nrank){
##  setTxtProgressBar(pb, jj, label=paste( round(jj/nrank*100, 0),"% done"))
##  
##    rf.out.rel<-suppressWarnings(randomForest(id_rel~delta*ngal*pval*ind*p_sw*color*p_lillie,data=trainset,importance=TRUE))
##    #rf.out.rel<-randomForest(id_rel~delta*ngal*pval*ind*p_sw*color*p_lillie*sph*sph.dens*gap.sph,data=trainset,importance=TRUE)
##    rf.pred.rel<-predict(rf.out.rel,newdata=dat_cum)
##    rf.pred.aux<-rf.pred.aux+rf.pred.rel
##  
##    if(length(dat_cum$rf.pred.rel) > 0){
##      dat_cum$rf.pred.rel=rf.pred.rel
##      dat<-dat_cum
##    } else {
##      dat<-data.frame(dat_cum,rf.pred.rel)
##    }
##    dat<-subset(dat,dat$rf.pred.rel>lim)
##   
##    if(length(dat$ngroup)>0){
##      cont.rel=cont.rel+1
##      if(cont.rel==1){
##        results=dat
##      } else {
##        results<-rbind(results,dat)
##      }
##    }
##  }
##  close(pb)
##  
##  if(cont.rel>0){
##    ntodo=length(dat_cum$ngroup)
##    rank<-1:2
##    for(i in 1:ntodo){
##      aux<-subset(results,results$ngroup == dat_cum$ngroup[i])
##      rank[i]=length(aux$ngroup)/nrank
##    }
##    rel_groups<-data.frame(dat_cum,rank)
##    if(length(dat_cum$ngroup)==1){
##      rel_groups<-rel_groups[1,]
##    }
##  } else {
##    rel_groups<-'No hay cumulos relajados'
##  }
##  
##  write.table(rel_groups,file=relaxed.name,row.names=FALSE)
##  rf.pred.rel<-rf.pred.aux/nrank
##  if(length(dat_cum$rf.pred.rel)>0){
##    dat_cum$rf.pred.rel=rf.pred.rel
##  } else {
##    dat_cum<-data.frame(dat_cum,rf.pred.rel)
##  }
##  write.table(dat_cum,file=name.groups,row.names=FALSE)
  #}}}

 
