#========================
#  strain dominance
# require(cgam)

require(Iso)
require(lubridate)

rm( list=ls() )
set.seed(523)

###########################
#  set-up parallel computations
#  uses forking
library(parallel)
MaxCores <- detectCores() -  1


########################################################################################

# source functions
source("myfunctions_v3.r")

# get data

all.files <- list.files()
data.files <- all.files[ grep( ".csv", all.files ) ]
dat.tmp <- read.table(data.files[1], header=T, sep=",", as.is = T)
dat.tmp[,2] <- as.Date(dat.tmp[,2])
dat.cnt     <- split( dat.tmp, dat.tmp[,1] )


#===============================================================================

#  number of bootstrap samples
N.perm <- 1000


#======================================================================================
# B.1.1.7

first.day.b117 <- as.Date("2020-09-20")
country.b117 <- match(c("Netherlands","United Kingdom"), names(dat.cnt) )

dat <- dat.cnt[country.b117]

lastday.b117 <- vector("list",length(dat))
pval.b117 <- lastday.b117
for ( k in 1:length(dat) ) {
  nl <- dim(dat[[k]])[1]
  sl <- sum(dat[[k]][,2] < first.day.b117 )
  lastday.b117[[k]] <- dat[[k]][(sl+7):nl,2]
  pval.b117[[k]] <- rep(-1,length(lastday.b117[[k]]))
}

for ( k in 1:length(dat) ){
  
  dd <- dat[[k]][,c(2,3,4)]
  dd[,2] <- dd[,2]-dd[,3]
  dd[,1] <- as.Date(dd[,1])
  
  for ( kk in 1:length(lastday.b117[[k]] ) ){
    print(paste(kk,k,names(dat)[k],sep="   "))
    # time range
  
    idx.window <- ( dd[,1] >= first.day.b117 ) & ( dd[,1] <= lastday.b117[[k]][kk])
    da <- dd[idx.window,]
    if ( sum(da[,3]) < 10 ) next 
    
    # baseline
    ffit.constant <- glm( as.matrix(da[,c(3,2)]) ~ 1, family = binomial ) # null model (constant)
    ffit.iso      <- IsoLog(da[,c(3,2)], decreasing = F )
    
    # calculate probabilities
    p2 <- da[,3]/(da[,2]+da[,3])         # full model
    p1 <- exp(ffit.iso)/(1+exp(ffit.iso))  # isotonic model
    p0 <- sum(da[,3])/sum(da[,2]+da[,3])            # null model
    
    # calculate loglikelihood ratio
    dloglik <- sum( da[,3]*log(p1/p0) + da[,2]*log((1-p1)/(1-p0)) )  # loglikelihood ratio
    
    
    #  expand for bootstrap
    ddl <- as.data.frame( t( da[,2:3] ) )
    DD <- cbind( rep( as.character(da[,1]), da[,2]+da[,3]),
                 unlist( lapply( ddl, FUN=function(x) rep(c("W","M"),x) ) )  )
    
    #  do the bootstrap
    zz <- vector("list", length=N.perm)
    LL <- unlist( mclapply(zz,
                           FUN=Boot.Pval,
                           DD[,2],
                           DD[,1],
                           p0, 
                           Decreasing=F,
                           Labels=c("M","W"),
                           mc.cores = MaxCores) )
    
    
    # calculate the p-value
    pval.b117[[k]][kk] <- ( sum( dloglik <= LL ) +1) /( N.perm + 2 )
  }
  
}



#======================================================================================
# B.1.351

first.day.b1351 <- as.Date("2020-10-01")
country.b1351 <- match(c("Netherlands"), names(dat.cnt) )

dat <- dat.cnt[country.b1351]

lastday.b1351 <- vector("list",length(dat))
pval.b1351 <- lastday.b1351
for ( k in 1:length(dat) ) {
  nl <- dim(dat[[k]])[1]
  sl <- sum(dat[[k]][,2] < first.day.b1351 )
  lastday.b1351[[k]] <- dat[[k]][(sl+7):nl,2]
  pval.b1351[[k]] <- rep(-1,length(lastday.b1351[[k]]))
}

for ( k in 1:length(dat) ){
  
  dd <- dat[[k]][,c(2,3,5)]
  dd[,2] <- dd[,2]-dd[,3]
  dd[,1] <- as.Date(dd[,1])
  
  for ( kk in 1:length(lastday.b1351[[k]] ) ){
    print(paste(kk,k,names(dat)[k],sep="   "))
    # time range
    
    idx.window <- ( dd[,1] >= first.day.b1351 ) & ( dd[,1] <= lastday.b1351[[k]][kk])
    da <- dd[idx.window,]
    if ( sum(da[,3]) < 10 ) next 
    
    # baseline
    ffit.constant <- glm( as.matrix(da[,c(3,2)]) ~ 1, family = binomial ) # null model (constant)
    ffit.iso      <- IsoLog(da[,c(3,2)], decreasing = F )
    
    # calculate probabilities
    p2 <- da[,3]/(da[,2]+da[,3])         # full model
    p1 <- exp(ffit.iso)/(1+exp(ffit.iso))  # isotonic model
    p0 <- sum(da[,3])/sum(da[,2]+da[,3])            # null model
    
    # calculate loglikelihood ratio
    dloglik <- sum( da[,3]*log(p1/p0) + da[,2]*log((1-p1)/(1-p0)) )  # loglikelihood ratio
    
    
    #  expand for bootstrap
    ddl <- as.data.frame( t( da[,2:3] ) )
    DD <- cbind( rep( as.character(da[,1]), da[,2]+da[,3]),
                 unlist( lapply( ddl, FUN=function(x) rep(c("W","M"),x) ) )  )
    
    #  do the bootstrap
    zz <- vector("list", length=N.perm)
    LL <- unlist( mclapply(zz,
                           FUN=Boot.Pval,
                           DD[,2],
                           DD[,1],
                           p0, 
                           Decreasing=F,
                           Labels=c("M","W"),
                           mc.cores = MaxCores) )
    
    
    # calculate the p-value
    pval.b1351[[k]][kk] <- ( sum( dloglik <= LL ) +1) /( N.perm + 2 )
  }
  
}



#======================================================================================
# B.1.1.7

first.day.b117 <- as.Date("2020-09-20")
country.b117 <- match(c("Netherlands","United Kingdom"), names(dat.cnt) )

dat <- dat.cnt[country.b117]

lastday.b117 <- vector("list",length(dat))
pval.b117 <- lastday.b117
for ( k in 1:length(dat) ) {
  nl <- dim(dat[[k]])[1]
  sl <- sum(dat[[k]][,2] < first.day.b117 )
  lastday.b117[[k]] <- dat[[k]][(sl+7):nl,2]
  pval.b117[[k]] <- rep(-1,length(lastday.b117[[k]]))
}

for ( k in 1:length(dat) ){
  
  dd <- dat[[k]][,c(2,3,4)]
  dd[,2] <- dd[,2]-dd[,3]
  dd[,1] <- as.Date(dd[,1])
  
  for ( kk in 1:length(lastday.b117[[k]] ) ){
    print(paste(kk,k,names(dat)[k],sep="   "))
    # time range
    
    idx.window <- ( dd[,1] >= first.day.b117 ) & ( dd[,1] <= lastday.b117[[k]][kk])
    da <- dd[idx.window,]
    if ( sum(da[,3]) < 10 ) next 
    
    # baseline
    ffit.constant <- glm( as.matrix(da[,c(3,2)]) ~ 1, family = binomial ) # null model (constant)
    ffit.iso      <- IsoLog(da[,c(3,2)], decreasing = F )
    
    # calculate probabilities
    p2 <- da[,3]/(da[,2]+da[,3])         # full model
    p1 <- exp(ffit.iso)/(1+exp(ffit.iso))  # isotonic model
    p0 <- sum(da[,3])/sum(da[,2]+da[,3])            # null model
    
    # calculate loglikelihood ratio
    dloglik <- sum( da[,3]*log(p1/p0) + da[,2]*log((1-p1)/(1-p0)) )  # loglikelihood ratio
    
    
    #  expand for bootstrap
    ddl <- as.data.frame( t( da[,2:3] ) )
    DD <- cbind( rep( as.character(da[,1]), da[,2]+da[,3]),
                 unlist( lapply( ddl, FUN=function(x) rep(c("W","M"),x) ) )  )
    
    #  do the bootstrap
    zz <- vector("list", length=N.perm)
    LL <- unlist( mclapply(zz,
                           FUN=Boot.Pval,
                           DD[,2],
                           DD[,1],
                           p0, 
                           Decreasing=F,
                           Labels=c("M","W"),
                           mc.cores = MaxCores) )
    
    
    # calculate the p-value
    pval.b117[[k]][kk] <- ( sum( dloglik <= LL ) +1) /( N.perm + 2 )
  }
  
}



#======================================================================================
# R.1

first.day.r1 <- as.Date("2020-10-24")
last.day.r1 <- as.Date(c("2021-02-15",
            "2021-02-28",
            "2021-03-15",
            "2021-03-25",
            "2021-04-10" ) )

dat <- dat.cnt
pval.r1 <- matrix(0,length(dat),length(last.day.r1))
rownames(pval.r1) <- names(dat)
colnames(pval.r1) <- as.character( last.day.r1 )

for ( kk in 1:length(last.day.r1) ){
  for ( k in 1:length(dat) ){
  
   # select country and strain
    print(paste(kk,k,names(dat)[k],sep="   "))
    dd <- dat[[k]][,c(2,3,6)]
    dd[,2] <- dd[,2]-dd[,3]
  
    # time range
    idx.window <- ( dd[,1] >= first.day.r1 ) & ( dd[,1] <= last.day.r1[[kk]])
    print( paste(names(dat)[k],sum(idx.window),sep="    " ) )
    pval.r1[k,kk] <- -2
    if ( sum(idx.window ) <=3 ) next
    dd <- dd[idx.window,]
  
    # exclusion filter 
    ss1 <- sum(dd[,3])
    ss2 <- sum(dd[,2])
    pval.r1[k,kk] <- -1
    print( paste( names(dat)[k], ss1, ss2, sep="    ") )
    if (( ss1 < 20 )||(ss2 < 20 )) next 
  
    # baseline
    ffit.constant <- glm( as.matrix(dd[,c(3,2)]) ~ 1, family = binomial ) # null model (constant)
    ffit.iso      <- IsoLog(dd[,c(3,2)], decreasing = F )
  
    # calculate probabilities
    p2 <- dd[,3]/(dd[,2]+dd[,3])                    # full model
    p1 <- exp(ffit.iso)/(1+exp(ffit.iso))           # isotonic model
    p0 <- sum(dd[,3])/sum(dd[,2]+dd[,3])            # null model
  
    # calculate loglikelihood ratio
    dloglik <- sum( dd[,3]*log(p1/p0) + dd[,2]*log((1-p1)/(1-p0)) )  # loglikelihood ratio
  
  
    #  expand for bootstrap
    ddl <- as.data.frame( t( dd[,2:3] ) )
    DD <- cbind( rep( as.character(dd[,1]), dd[,2]+dd[,3]),
                unlist( lapply( ddl, FUN=function(x) rep(c("W","M"),x) ) )  )
  
    #  do the bootstrap
    zz <- vector("list", length=N.perm)
    LL <- unlist( mclapply(zz,
                           FUN=Boot.Pval,
                           DD[,2],
                           DD[,1],
                           p0, 
                           Decreasing=F,
                           Labels=c("M","W"),
                           mc.cores = MaxCores) )
  
  
    # calculate the p-value
    pval.r1[k,kk] <- ( sum( dloglik <= LL ) +1) /( N.perm + 2 )
  
  }
}

idx <- apply(pval.r1 < -0.5,1,sum ) < length(last.day.r1)
Pval.r1 <- pval.r1[idx,]


#=============================================================================================
#
dat.tmp <- read.table(data.files[2],sep=",",header = T, as.is = T)

dat.cnt     <- split( dat.tmp[,1:4], dat.tmp[,1] )

FF.include  <- function(x){
  c1 <- dim(x)[1]  > 9             # at least 10 days
  s1 <- sum( x[,4] )
  s2 <- sum( x[,3] )
  c2 <- (s2 - s1) > 19             # at least 20 cases of other
  c3 <- s1 > 19                    # at least 20 cases of variant
  return( c1 & c2 & c3 )
}

idx.d614g   <- sapply( dat.cnt, FF.include )
dat.cnt.d614g <- dat.cnt[ idx.d614g ]

#======================================================================================
# D614G

dat <- dat.cnt.d614g
pval.d614g <- rep(0,length(dat))
names(pval.d614g) <- names(dat)

for ( k in 1:length(dat) ){
  print(paste(k,names(dat)[k],sep="   "))
  dd     <- dat[[k]][,c(2,3,4)]
  iidx   <- dd[,2] > 0
  dd     <- dd[ iidx, ]
  dd[,2] <- dd[,2]-dd[,3]

  # baseline
  ffit.constant <- glm( as.matrix(dd[,c(3,2)]) ~ 1, family = binomial ) # null model (constant)
  ffit.iso      <- IsoLog( as.matrix(dd[,c(3,2)]), decreasing = F )
  
  # calculate probabilities
  p2 <- dd[,3]/(dd[,2]+dd[,3])                    # full model
  p1 <- exp(ffit.iso)/(1+exp(ffit.iso))           # isotonic model
  p0 <- sum(dd[,3])/sum(dd[,2]+dd[,3])            # null model
  
  # calculate loglikelihood ratio
  dloglik <- sum( dd[,3]*log(p1/p0) + dd[,2]*log((1-p1)/(1-p0)) )  # loglikelihood ratio
  
  
  #  expand for bootstrap
  ddl <- as.data.frame( t( dd[,2:3] ) )
  DD <- cbind( rep( as.character(dd[,1]), dd[,2]+dd[,3]),
               unlist( lapply( ddl, FUN=function(x) rep(c("W","M"),x) ) )  )
  
  #  do the bootstrap
  zz <- vector("list", length=N.perm)
  LL <- unlist( mclapply(zz,
                         FUN=Boot.Pval,
                         DD[,2],
                         DD[,1],
                         p0, 
                         Decreasing=F,
                         Labels=c("M","W"),
                         mc.cores = MaxCores) )
  
  
  # calculate the p-value
  pval.d614g[k] <- ( sum( dloglik <= LL ) +1) /( N.perm + 2 )
  
}


