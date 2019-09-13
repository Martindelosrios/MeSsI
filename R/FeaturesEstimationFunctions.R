# DresslerShectmanTest
#{{{
#' DresslerShectmanTest
#' @description This function returns the value of the Dressler-Shectman statistic for a galaxy cluster.
#' @param 
#' cat data frame with a the angular coordinates (ra and dec in radians) and the line of sight velocity (vel in km/s). The columns must be named 'ra', 'dec' and 'vel'
#' @return The Delta value of the Dressler-Shectman test corresponding to the galaxy cluster.
#' @export
#' @examples
#' DresslerShectmanTest(cat)

DresslerShectmanTest <- function(cat){
  n <- length(cat[,1]) # Number of galaxies in the catalog
  
  x   <- cat$ra # x coordinate
  y   <- cat$dec # y coordinate
  vel <- cat$vel # velocity

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
#' cat data frame with a the angular coordinates (ra and dec in radians) and the line of sight velocity (vel in km/s). The columns must be named 'ra', 'dec' and 'vel'
#' @return The Delta value of the Dressler-Shectman test corresponding to the galaxy cluster.
#' @export
#' @examples
#' DresslerShectmanTest2(cat)

DresslerShectmanTest2 <- function(cat){
  n <- length(x) # Number of galaxies in the catalog
  
  x   <- cat$ra # x coordinate 
  y   <- cat$dec # y coordinate
  vel <- cat$vel # velocity

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
#' cat data frame with a the angular coordinates (ra and dec in radians) and the line of sight velocity (vel in km/s). The columns must be named 'ra', 'dec' and 'vel'
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
            vel      <- cat$vel
            r        <- get_distance(ra, x0, dec, y0, ngal = floor(sqrt(ra)))
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
#' cat data frame with a the angular coordinates (ra and dec in radians) and the line of sight velocity (vel in km/s). The columns must be named 'ra', 'dec' and 'vel'
#' @export
#' @examples
#' DresslerShectmanPval(cat)

DresslerShectmanPval <- function(cat){
  n     <- length(cat$vel)
  nran  <- 0
  nmont <- 100
  del   <- DresslerShectmanTest(cat)
  mat   <- cat # we will modify mat in each iteration

  for(i in 1:nmont){
    mat$vel <- sample(mat$vel, size = n, replace = FALSE)
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
#' cat data frame with a the angular coordinates (ra and dec in radians) and the line of sight velocity (vel in km/s). The columns must be named 'ra', 'dec' and 'vel'
#' @export
#' @examples
#' DresslerShectmanPval2(cat)

DresslerShectmanPval2 <- function(cat){
  n     <- length(cat$vel)
  nran  <- 0
  nmont <- 100
  del   <- DresslerShectmanTest(cat)
  mat   <- cat # we will modify mat in each iteration

  for(i in 1:nmont){
    mat$vel <- sample(mat$vel, size = n, replace = FALSE)
    del_mon <- DresslerShectmanTest2(mat)
    if(del_mon > del){
      nran <- nran+1
    }
  }
  nran <- nran/nmont
  return(nran)
}
#}}}
#DresslerShectmanIndividual 
#{{{
#' DresslerShectmanIndividual
#' @description This function returns the delta value of the Dressler-Shectman statistic each galaxy of a galaxy cluster with ngal = 10.
#' @return The delta value of the Dressler-Shectman test corresponding to each galaxy of a galaxy cluster.
#' @param 
#' group data frame with a the angular coordinates (ra and dec in radians) and the line of sight velocity (vel in km/s). The columns must be named 'ra', 'dec' and 'vel'.
#' @export
#' @examples
#' DresslerShectmanIndividual(group)

DresslerShectmanIndividual <- function(group){
  n <- length(group$ra)
  # Estimation of the delta statistic for each galaxy 
  delta <- 1:n
  for(j in 1:n){
    r        <- get_distance(dat$ra, dat$ra[j], dat$dec, dat$dec[j], ngal  = 10)
    delta[j] <- DresslerShectmanGalaxy(dat$vel, r)
  }

  return(delta)
}
#}}}
#DresslerShectmanIndividual2 
#{{{
#' DresslerShectmanIndividual2
#' @description This function returns the delta value of the Dressler-Shectman statistic each galaxy of a galaxy cluster with ngal = floor(sqrt(ngal)).
#' @return The delta value of the Dressler-Shectman test corresponding to each galaxy of a galaxy cluster.
#' @param 
#' group data frame with a the angular coordinates (ra and dec in radians) and the line of sight velocity (vel in km/s). The columns must be named 'ra', 'dec' and 'vel'.
#' @export
#' @examples
#' DresslerShectmanIndividual2(group)

DresslerShectmanIndividual2 <- function(group){
  n <- length(group$ra)
  # Estimation of the delta statistic for each galaxy 
  delta <- 1:n
  for(j in 1:n){
    r        <- get_distance(dat$ra, dat$ra[j], dat$dec, dat$dec[j], ngal  = floor(sqrt(n)))
    delta[j] <- DresslerShectmanGalaxy(dat$vel, r)
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
#' @description This function estimates the projected distances between galaxies needed for the estiamtion of the virial radius of the galaxy cluster.
#' @param 
#' cat data frame with a the angular coordinates (ra and dec in radians) and the line of sight velocity (vel in km/s). The columns must be named 'ra', 'dec' and 'vel'
#' @return sum(1/d) where d is a vector with the projected distances.
#' @export
#' @examples
#' ProjectedDistance(cat)

ProjectedDistance <- function(cat){
  id  <- sort(cat$vel, decreasing = FALSE, index.return=TRUE)$ix
  vel <- cat$vel[id]
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
            r[d] <- r[d]*D.A(mean(vel/300000))
          }
        }
      }
    }
  }
  r <- 1/r
  r <- subset(r,r<10000000)
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
  ngal.lim <- 30
  if(length(ra) > ngal.lim){ # Cut in the number of member galaxies
    ngroup       <- dat$id[1]
    ra           <- mean(dat$ra)
    dec          <- mean(dat$dec)
    z            <- mean(dat$z)
    ngal         <- length(dat$ra)
    color        <- mean(dat$color)
    sort_mags    <- sort(group$mag, decreasing = FALSE)
    mag_max      <- sort_mags[1]
    gap_mag      <- sort_mags[2] - sort_mags[1]
    Delta        <- DresslerShectmanTest(dat$ra, dat$dec, dat$z*c)
    Delta2       <- DresslerShectmanTest2(dat$ra, dat$dec, dat$z*c)
    pval_ds      <- DresslerShectmanPval(dat$ra, dat$dec, dat$z*c)
    ind          <- DresslerShectmanIter(dat$ra, dat$dec, dat$z*c)    
    pval_sw      <- shapiro.test(dat$z*c)$p.val
    pval_sf      <- sf.test(dat$z*c)$p.val
    pval_ad      <- ad.test(dat$z*c)$p.val
    pval_cvm     <- suppressWarnings(cvm.test(vel))$p.val
    pval_lillie  <- lillie.test(dat$z*c)$p.val
    pval_pearson <- pearson.test(dat$z*c)$p.val
  }  
  features <- data.frame(ngroup, ra, dec, z, ngal, color, , mag_max, gap_mag, Delta, Delta2, pval_ds, ind, pval_sw, pval_sf, pval_ad, pval_cvm, pval_lillie, pval_pearson)
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

get_cluster_features <- function(dat){
  for(i in 1:ntotal){
    if(length(dat$ra)>0){
      groupid  <- dat$id[1] # Id if the group that will be studied
      group    <- subset(dat,dat$id==groupid) # Select the galaxies of the group
      dat      <- subset(dat,dat$id != groupid) # Remove the galaxies of the group from the general catalog
      features <- ClusterFeatures(group, ...) # Estimates the features
      if(exists('Allfeatures') == FALSE){Allfeatures <- features}
      features <- rbind(Allfeatures, features)
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
  vel <- group_gals$z*c
  
  delta  <- DresslerShectmanIndividual(x, y, vel)
  delta2 <- DresslerShectmanIndividual2(x, y, vel)
  
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
    pval_lillie[j]  <- lillie.test(vel_loc)p.value
    pval_pearson[j] <- pearson.test(vel_loc)$p.value
  }
  
  ngroup      <- group_data$ngroup
  ra          <- x
  dec         <- y
  z           <- vel/c
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

get_galaxies_features(dat, ClustersData){
  ncluster <- length(ClustersData$ngroup) # Total number of galaxy clusters

  for(i in 1:nclusters){
    group_gals    <- subset(dat, dat$id == ClustersData$ngroup[i]) # Id of the cluster that will be studied in this iteration
    group_data    <- ClustersData[i,] # Data of the cluster that will be studied in this iteration
    features <- GalaxiesFeatures(group_gals, group_data) # Estimation of the galaxy features
    if(exists('Allfeatures') == FALSE){Allfeatures <- features}
    features <- rbind(Allfeatures, features)
  }
  write.table(Allfeatures, file = name.gal, row.names = FALSE, quote = FALSE)
}

#}}}




#messi
#{{{
#dat debe ser un data frame con: id (numero),ra[° decimales],dec[° decimales],z (redshift),mag(aparente),color

messi <- function(dat, cum = TRUE, gal = TRUE, rank = TRUE, relaxed = TRUE, 
                  name.groups = 'estadisticos_grupos.dat',
                  name.gal = 'estadisticos_galaxias.dat',
                  rank.name = 'ranking.dat', relaxed.name = 'ranking_relaxed.dat',
                  folder = 'folder', nrank = 100, name_trainset_cum = 'trainset_cum.dat',
                  name_trainset_gal = 'trainset_gal.dat', ntotal = 0, 
                  ngal.lim = 30, est = TRUE, ...){

  # Configuration and storage folders
#{{{
  if(file.exists(folder)==FALSE){
    system(paste('mkdir',folder,sep=' '))
  }
  mc <- paste(folder,'/merging_clusters',sep='') 
  if(file.exists(mc)==FALSE){
    system(paste('mkdir',mc,sep=' '))
  }
  name.groups  <- paste(folder,'/',name.groups,sep='')
  name.gal     <- paste(folder,'/',name.gal,sep='')
  rank.name    <- paste(folder,'/',rank.name,sep='')
  relaxed.name <- paste(folder,'/',relaxed.name,sep='')
#}}}  
  t1 <- proc.time()

  if(cum == TRUE){
    print('Starting the estimation of the galaxy clusters features')
    if(ntotal==0){ntotal <- length(dat$ra)}
    ClustersData <- get_cluster_features(dat)
  }
  if(gal == TRUE){
    print('Starting the estimation of the galaxy features')
    GalaxiesData <- get_galaxies_features(dat, ClustersData)
  }
  if(rank == TRUE){
    print('Starting the classification of the galaxy clusters')
    ClustersClassification <- get_clusters_classification(ClustersData)
  }
  if(est == TRUE){
    estabilidad_cum(folder = folder, n_it = 10, par = FALSE)
  }

  t2 <- proc.time()
  t  <- (t2[3]-t1[3])
  print(paste('The program delay',toString(t/60),'minutes',sep=' '))
}
#}}}

    # RANKING-RELAXED
  #{{{
  trainset<-read.table(file=name_trainset_cum,header=TRUE)
  nrank=nrank
  
  dat_cum<-read.table(file=name.groups,header=TRUE)
  if(length(dat_cum$rf.pred.rel) > 0){ #Si ya existe la medicion aca la elimino para luego reemplazarla
    dat_cum <- subset(dat_cum, select = -c(rf.pred.rel))
  }
  lim=0.6
  
  # Lectura de datos de cumulos
  #if(mean(dat_cum$z) > 0.15){
  #  trainset_rel<-read.table('/media/martin/store1/trabajos/mock/guo/mock_cumulos_snap56.dat',header=TRUE)
  #} else {
  #  trainset_rel<-read.table('/media/martin/store1/trabajos/mock/guo/mock_cumulos_snap63.dat',header=TRUE)
  #}
  cont.rel=0
  rf.pred.aux<-1:length(dat_cum$ngroup)
  rf.pred.aux[]=0
  
  pb <- txtProgressBar(title = "progress bar", min = 0,max = nrank, width = 82)
  for(jj in 1:nrank){
  setTxtProgressBar(pb, jj, label=paste( round(jj/nrank*100, 0),"% done"))
  
    rf.out.rel<-suppressWarnings(randomForest(id_rel~delta*ngal*pval*ind*p_sw*color*p_lillie,data=trainset,importance=TRUE))
    #rf.out.rel<-randomForest(id_rel~delta*ngal*pval*ind*p_sw*color*p_lillie*sph*sph.dens*gap.sph,data=trainset,importance=TRUE)
    rf.pred.rel<-predict(rf.out.rel,newdata=dat_cum)
    rf.pred.aux<-rf.pred.aux+rf.pred.rel
  
    if(length(dat_cum$rf.pred.rel) > 0){
      dat_cum$rf.pred.rel=rf.pred.rel
      dat<-dat_cum
    } else {
      dat<-data.frame(dat_cum,rf.pred.rel)
    }
    dat<-subset(dat,dat$rf.pred.rel>lim)
   
    if(length(dat$ngroup)>0){
      cont.rel=cont.rel+1
      if(cont.rel==1){
        results=dat
      } else {
        results<-rbind(results,dat)
      }
    }
  }
  close(pb)
  
  if(cont.rel>0){
    ntodo=length(dat_cum$ngroup)
    rank<-1:2
    for(i in 1:ntodo){
      aux<-subset(results,results$ngroup == dat_cum$ngroup[i])
      rank[i]=length(aux$ngroup)/nrank
    }
    rel_groups<-data.frame(dat_cum,rank)
    if(length(dat_cum$ngroup)==1){
      rel_groups<-rel_groups[1,]
    }
  } else {
    rel_groups<-'No hay cumulos relajados'
  }
  
  write.table(rel_groups,file=relaxed.name,row.names=FALSE)
  rf.pred.rel<-rf.pred.aux/nrank
  if(length(dat_cum$rf.pred.rel)>0){
    dat_cum$rf.pred.rel=rf.pred.rel
  } else {
    dat_cum<-data.frame(dat_cum,rf.pred.rel)
  }
  write.table(dat_cum,file=name.groups,row.names=FALSE)
  #}}}
    # RANKING
  #{{{
  nrank=nrank
  
  lim=0.3
  lim.gal=0.4
  #lim=0.2 #Version2
  #lim.gal=0.3 #Version2
  
  trainset<-read.table(file=name_trainset_cum,header=TRUE)
  trainset_gal<-read.table(file=name_trainset_gal,header=TRUE)
  
  
  dat_gals<-read.table(file=name.gal,header=TRUE)
  if(length(dat_gals$rf.pred.gal) > 0){ #Si ya existe la medicion aca la elimino para luego reemplazarla
    dat_gals <- subset(dat_gals, select = -c(rf.pred.gal))
  }
  dat_cum<-read.table(file=name.groups,header=TRUE)
  if(length(dat_cum$rf.pred.mer) > 0){ #Si ya existe la medicion aca la elimino para luego reemplazarla
    dat_cum <- subset(dat_cum, select = -c(rf.pred.mer))
  }
  
  cont.results=0
  rf.pred.aux<-1:length(dat_cum$ngroup)
  rf.pred.aux[]=0
  
  pb <- txtProgressBar(title = "progress bar", min = 0,max = nrank, width = 82)
  for(jj in 1:nrank){
  setTxtProgressBar(pb, jj, label=paste( round(jj/nrank*100, 0),"% done"))
    
    rf.out<-suppressWarnings(randomForest(id_mer~delta*ngal*pval*ind*p_sw*color*p_lillie,data=trainset,importance=TRUE))
    rf.pred.mer<-predict(rf.out,newdata=dat_cum)
    rf.pred.aux<-rf.pred.aux+rf.pred.mer
    
    if(length(dat_cum$rf.pred.mer) > 0){
      dat_cum$rf.pred.mer=rf.pred.mer
      dat<-dat_cum
    } else {
      dat<-data.frame(dat_cum,rf.pred.mer)
    }
    dat<-subset(dat,dat$rf.pred.mer>lim)
    
    
    rf.out.gal<-suppressWarnings(randomForest(id~deltas_gal*del_cum*g_r*sw_gal*col_cum,data=trainset_gal,importance=TRUE))
    rf.pred.gal<-predict(rf.out.gal,newdata=dat_gals)
    
    if(length(dat_gals$rf.pred.gal) > 0){
      dat_gals$rf.pred.gal=rf.pred.gal
      dat_gal<-dat_gals
    } else {
      dat_gal<-data.frame(dat_gals,rf.pred.gal)
    }
    dat_gal<-subset(dat_gal,dat_gal$rf.pred.gal>lim.gal)
    
    
    ncum=length(dat$delta)
    
    
    if(ncum > 0){ 
      aux<-1:(ncum*22)
      aux[1:(ncum*22)]=-99
      mat<-matrix(aux,ncol=22,nrow=ncum)
      for(i in 1:ncum){
        group<-subset(dat_gal,dat_gal$ngroup==dat$ngroup[i])
        if(length(group$ra)>10){
       
        ra=group$ra
        dec=group$dec
        vel=group$redshift*300000
        mag=group$mag_r
        color=group$g_r
        delta=group$deltas_gal
        peso=group$rf.pred.gal 
        
        delmin=0
        grupo.id<-paste(folder,'/merging_clusters/',toString(group$ngroup[1]),sep='')
        mixt.real(ra,dec,vel,mag,color,delta,peso,delmin,grupo.id)->v
        mat[i,1]=dat$ngroup[i]
        mat[i,2:22]=v
        }
      }    
      mat<-data.frame(mat)
      vec<-c('ngroup','l1','l2','rvir1','rvir2','dvel1','dvel2','mas1','mas2','par1','par2','tot','sigma1_ra','sigma1_dec','sigma2_ra','sigma2_dec','racen1','racen2','deccen1','deccen2','velcen1','velcen2')
      colnames(mat)<-vec
      mat<-subset(mat,mat$mas1>0)
      #mat<-subset(mat,mat$par1>0.22)
      #mat<-subset(mat,mat$par1>0.15)
     
      cont.results=cont.results+1
      if(cont.results==1){
        results=mat$ngroup
        results_mat=mat
      } else {
        results<-c(results,mat$ngroup)
        results_mat<-rbind(results_mat,mat)
      }
    } 
  }
  close(pb)
  
  
  if(cont.results > 0){
    ntodo=length(results)
    name<-1:2
    largo<-1:2
    c=0
    for(i in 1:ntodo){
       if(length(results)>0){
          s1<-subset(results,results == results[1])
          c=c+1
          name[c]=s1[1]
          largo[c]=length(s1)
          results<-subset(results,results != results[1])
       }
    }
    
    largo=largo/nrank
    ngroup<-1:length(largo)
    l1<-1:length(largo)
    l2<-1:length(largo)
    l1_err<-1:length(largo)
    l2_err<-1:length(largo)
    rvir1<-1:length(largo)
    rvir2<-1:length(largo)
    rvir1_err<-1:length(largo)
    rvir2_err<-1:length(largo)
    dvel1<-1:length(largo)
    dvel2<-1:length(largo)
    dvel1_err<-1:length(largo)
    dvel2_err<-1:length(largo)
    m1<-1:length(largo)
    m2<-1:length(largo)
    par1<-1:length(largo)
    par2<-1:length(largo)
    par1_err<-1:length(largo)
    par2_err<-1:length(largo)
    tot<-1:length(largo)
    tot_err<-1:length(largo)
    sigma1_ra<-1:length(largo) 
    sigma2_ra<-1:length(largo) 
    sigma1_ra_err<-1:length(largo) 
    sigma2_ra_err<-1:length(largo) 
    sigma1_dec<-1:length(largo) 
    sigma2_dec<-1:length(largo) 
    sigma1_dec_err<-1:length(largo) 
    sigma2_dec_err<-1:length(largo) 
    m1_err<-1:length(largo)
    m2_err<-1:length(largo)
    ra1<-1:length(largo)
    ra2<-1:length(largo)
    ra1_err<-1:length(largo)
    ra2_err<-1:length(largo)
    dec1<-1:length(largo)
    dec2<-1:length(largo)
    dec1_err<-1:length(largo)
    dec2_err<-1:length(largo)
    z1<-1:length(largo)
    z2<-1:length(largo)
    z1_err<-1:length(largo)
    z2_err<-1:length(largo)
    
    for(i in 1:length(largo)){
      if(length(results_mat$ngroup)>0){
        mat_aux<-subset(results_mat,results_mat$ngroup==results_mat$ngroup[1])
        results_mat<-subset(results_mat,results_mat$ngroup!=results_mat$ngroup[1])
      
        ngroup[i]=mat_aux$ngroup[1]
        l1[i]=mean(mat_aux$l1) 
        l1_err[i]=sd(mat_aux$l1) 
        l2[i]=mean(mat_aux$l2) 
        l2_err[i]=sd(mat_aux$l2) 
        rvir1[i]=mean(mat_aux$rvir1) 
        rvir1_err[i]=sd(mat_aux$rvir1) 
        rvir2[i]=mean(mat_aux$rvir2) 
        rvir2_err[i]=sd(mat_aux$rvir2) 
        dvel1[i]=mean(mat_aux$dvel1) 
        dvel1_err[i]=sd(mat_aux$dvel1) 
        dvel2[i]=mean(mat_aux$dvel2) 
        dvel2_err[i]=sd(mat_aux$dvel2) 
        m1[i]=mean(mat_aux$mas1) 
        m1_err[i]=sd(mat_aux$mas1)
        m2[i]=mean(mat_aux$mas2) 
        m2_err[i]=sd(mat_aux$mas2)
        par1[i]=mean(mat_aux$par1) 
        par1_err[i]=sd(mat_aux$par1)
        par2[i]=mean(mat_aux$par2) 
        par2_err[i]=sd(mat_aux$par2)
        tot[i]=mean(mat_aux$tot)
        sigma1_ra[i]=mean(mat_aux$sigma1_ra) 
        sigma1_ra_err[i]=sd(mat_aux$sigma1_ra)
        sigma2_ra[i]=mean(mat_aux$sigma2_ra) 
        sigma2_ra_err[i]=sd(mat_aux$sigma2_ra)
        sigma1_dec[i]=mean(mat_aux$sigma1_dec) 
        sigma1_dec_err[i]=sd(mat_aux$sigma1_dec)
        sigma2_dec[i]=mean(mat_aux$sigma2_dec) 
        sigma2_dec_err[i]=sd(mat_aux$sigma2_dec)
        tot_err[i]=sd(mat_aux$tot)
        ra1[i]=mean(mat_aux$racen1)*180/pi
        ra1_err[i]=sd(mat_aux$racen1)*180/pi
        ra2[i]=mean(mat_aux$racen2)*180/pi
        ra2_err[i]=sd(mat_aux$racen2)*180/pi
        dec1[i]=mean(mat_aux$deccen1)*180/pi
        dec1_err[i]=sd(mat_aux$deccen1)*180/pi
        dec2[i]=mean(mat_aux$deccen2)*180/pi
        dec2_err[i]=sd(mat_aux$deccen2)*180/pi
        z1[i]=mean(mat_aux$velcen1)/300000
        z1_err[i]=sd(mat_aux$velcen1)/300000
        z2[i]=mean(mat_aux$velcen2)/300000
        z2_err[i]=sd(mat_aux$velcen2)/300000
      }
    }
    
    rank=largo
    rankin<-data.frame(ngroup,rank,m1,m1_err,ra1,ra1_err,dec1,dec1_err,z1,z1_err,m2,m2_err,ra2,ra2_err,dec2,dec2_err,z2,z2_err,l1,l1_err,l2,l2_err,tot,tot_err,rvir1,rvir1_err,rvir2,rvir2_err,dvel1,dvel1_err,dvel2,dvel2_err,par1,par1_err,par2,par2_err,sigma1_ra,sigma1_ra_err,sigma2_ra,sigma2_ra_err,sigma1_dec,sigma1_dec_err,sigma2_dec,sigma2_dec_err)
  
    if(c==1){
       rankin<-rankin[1,]
    }
  } else {
    rankin<-'No hay cumulos en merger'
  }
  write.table(rankin,file=rank.name,row.names=FALSE)
  rf.pred.mer<-rf.pred.aux/nrank
  if(length(dat_cum$rf.pred.mer)>0){
    dat_cum$rf.pred.mer=rf.pred.mer
  } else {
    dat_cum<-data.frame(dat_cum,rf.pred.mer)
  }
  write.table(dat_cum,file=name.groups,row.names=FALSE)
  #}}}
