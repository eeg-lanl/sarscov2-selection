#==============================================
#  functions
#----------------------------------------------
#  reworked function to speed up the bootstrap analysis


#----------------------------------------------------
#  Make Tables with specified categoriges
#  Enables fixed format output and zero entries
#----------------------------------------------------

Table <- function(x,Labels){
  #=c("D614","G614") ){
  xx <- rep(0, length(Labels) )
  names(xx) <- Labels
  xx.t <- table( x )
  xx[ names(xx.t) ] <- xx.t
  return(xx)
}

#=================================================================
#  Compute isotonic logistic regression, using the pava rountine
#  from Iso 
#  Input:  yy     matrix of successes, failures, ordered in time
#          decreasing = F
#          theta  initial value, default is zero
#          err    L1 error on estimated probs.  Recommended 10e-5
#          eps    numerical stability for prob going to 0 or 1
IsoLog <- function(Y,decreasing=F,theta=0,err=1e-5,min.w=1e-6){
  G <- 1
  theta0 <- theta
  pp <- exp(theta0)/(1+exp(theta0))
  
  yy   <- Y[,1]
  nn   <- apply(Y,1,sum)
  
  while ( G > err ){
    ww <- nn*pmax( pp*(1-pp), min.w )  # ensure the weights are always larger than eps
    zz <- theta0 + (yy-nn*pp)/ww
    theta1 <- pava(zz,ww,decreasing=decreasing)
    p1 <- exp(theta1)/(1+exp(theta1))
    G <- sum( abs(p1-pp) )
    #
    theta0 <- theta1
    pp <- p1
  }
  return(theta1)
}


#===================================================================================
#  fit a Ashaped 
Ashape <- function(Y,lmode,theta=0,err=1e-5,min.w=1e-6){
  G <- 1
  theta0 <- theta
  pp <- exp( theta0)/( 1+exp(theta0) )
  
  yy <- Y[,1]
  nn <- apply(Y,1,sum)
  iidx <- (1:length(yy)) < (lmode+1)

  while ( G > err ){
    ww <- nn*pmax( pp*(1-pp), min.w )  # ensure the weights are always larger than eps
    zz <- theta0 + (yy-nn*pp)/ww
    theta1 <- pava(zz[iidx],ww[iidx],decreasing=F)
    theta2 <- pava(zz[!iidx],ww[!iidx],decreasing=T)
    ttheta <- c(theta1,theta2)
    p1 <- exp(ttheta)/(1+exp(ttheta))
    G <- sum( abs(p1-pp) )
    #
    theta0 <- ttheta
    pp <- p1
  }
  loglik <- sum(yy*log(pp) + (nn-yy)*log(1-pp))
  return( list( theta=theta0,p=pp,loglik=loglik ) )
}

Ashape2 <- function(Y,Lmode=c(1,dim(Y)[1]),theta=0,err=1e-6,min.w=1e-6){
  lmode <- Lmode[1]:Lmode[2]
  Ffit <- vector("list",length = length(lmode))  
  Ffit[[1]] <- Ashape(Y,lmode[1],theta)
  for ( k in 2:length(lmode) ){
    Ffit[[k]] <- Ashape(Y,lmode=lmode[k],Ffit[[k-1]]$theta)
  }
  iidx <- which.max(sapply(Ffit, FUN=function(x) x$loglik))
  return( list( theta=Ffit[[iidx]]$theta, 
                p=Ffit[[iidx]]$p,
                loglik=Ffit[[iidx]]$loglik,
                mode = lmode[iidx] ) )
}

#===================================================================
#  Function to compute from one bootstrap sample
#  To be used in parallel
#  input: empty vector, needed to get it to run in parallel
#         y vector of 0-1
#         ddate: vector of dates
#         p0: probability under the null distribution
#         Decreasing: test decresing? T=yes, F=no
#         Labels: vector of labels to be used in Table
#===================================================================
Boot.Pval <- function(zz,yy,ddates,p0,Decreasing=F,Labels=c("M","W")){
  
    y.b <- sample(yy, replace = F)  # permutation of strain labels
    ddt.b <- split(y.b,ddates) # split by date
    BB.b <- t( sapply( ddt.b , Table, Labels=Labels ) )  
    
    ffit.b <- IsoLog( as.matrix( BB.b ), decreasing=Decreasing)
    p.b <- exp(ffit.b)/(1+exp(ffit.b))
    
    # calculate relative loglikelihood
    LL <- sum( BB.b[,1]*log(p.b/p0) + BB.b[,2]*log((1-p.b)/(1-p0)) )
    
    return(LL)
    
  }
