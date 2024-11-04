#########################################################################################
# This code is from 'SHT' r package source code with some modification for this paper.
#########################################################################################

rm(list=ls())
set.seed(1)
library(Rcpp)
library(MASS)
sourceCpp("cpp_mean2_mxPBF.cpp")
sourceCpp("cpp_cov2_mxPBF.cpp")
mean2.mxPBF <- function(X, Y, a0=0, b0=0, gamma=1.0, nthreads=1){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  if (ncol(X)!=ncol(Y)){
    stop("* mean2.mxPBF : two samples X and Y should be of same dimension.")
  }
  p = ncol(X)
  if ((nrow(X)<2)||(nrow(Y)<2)||(p<2)){
    stop("* mean2.mxPBF : inputs are invalid. Provide multivariate samples with multiple observations.")
  }
  if ((length(a0)>1)||(a0<0)){
    stop("* mean2.mxPBF : 'a0' should be a nonnegative number.")
  }
  if ((length(b0)>1)||(b0<0)){
    stop("* mean2.mxPBF : 'b0' should be a nonnegative number.")
  }
  if ((length(gamma)>1)||(gamma<=0)){
    stop("* mean2.mxPBF : 'gamma' should be a nonnegative number.")
  }
  nCores = as.integer(nthreads)
  if (nCores < 1){
    stop("* mean2.mxPBF : 'nthreads' should be a positive integer.")
  }
  
  ##############################################################
  # MAIN COMPUTATION
  if (nCores==1){
    log.BF.vec = as.vector(cpp_mean2_mxPBF_single(X, Y, a0, b0, gamma))  
  } else {
    log.BF.vec = as.vector(cpp_mean2_mxPBF_multiple(X, Y, a0, b0, gamma, nCores))
  }
  
  ##############################################################
  # FINALE
  hname   = "Two-sample Mean Test with Maximum Pairwise Bayes Factor"
  Ha      = "two means are not equal."
  
  thestat = max(log.BF.vec)
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "maximum BF"
  res   = list(statistic=thestat, alternative = Ha, method=hname, data.name = DNAME, log.BF.vec = log.BF.vec)
  class(res) = "htest"
  return(res)
}
cov2.mxPBF <- function(X, Y, a0=2.0, b0=2.0, gamma=1.0, nthreads=1){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  if (ncol(X)!=ncol(Y)){
    stop("* cov2.mxPBF : two samples X and Y should be of same dimension.")
  }
  p = ncol(X)
  if ((nrow(X)==1)||(nrow(Y)==1)||(p<2)){
    stop("* cov2.mxPBF : inputs are invalid. Provide multivariate samples with multiple observations.")
  }
  if ((length(a0)!=1)||(a0<=0)){
    stop("* cov2.mxPBF : 'a0' should be a nonnegative number.")
  }
  if ((length(b0)!=1)||(b0<=0)){
    stop("* cov2.mxPBF : 'b0' should be a nonnegative number.")
  }
  if ((length(gamma)!=1)||(gamma<=0)){
    stop("* cov2.mxPBF : 'gamma' should be a nonnegative number.")
  }
  nCores = as.integer(nthreads)
  if (nCores < 1){
    stop("* cov2.mxPBF : 'nthreads' should be a positive integer.")
  }
  
  ##############################################################
  # PRELIMINARY
  Xnew = as.matrix(scale(X, center=TRUE, scale=FALSE))
  Ynew = as.matrix(scale(Y, center=TRUE, scale=FALSE))
  
  ##############################################################
  # MAIN COMPUTATION
  if (nCores==1){
    log.BF.mat = cpp_cov2_mxPBF_single(Xnew, Ynew, a0, b0, gamma)
  } else {
    log.BF.mat = cpp_cov2_mxPBF_multiple(Xnew, Ynew, a0, b0, gamma, nCores)
  }
  diag(log.BF.mat) = -Inf
  
  ##############################################################
  # FINALE
  hname   = "Two-sample Covariance Test with Maximum Pairwise Bayes Factor"
  Ha      = "two covariances are not equal."
  
  thestat = max(log.BF.mat)
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "maximum BF"
  res   = list(statistic=thestat, alternative = Ha, method=hname, data.name = DNAME, log.BF.mat = log.BF.mat)
  class(res) = "htest"
  return(res)
}
check_nd <- function(x){
  cond1 = ((is.array(x))||(is.matrix(x)))
  cond2 = (all(!is.infinite(x)))
  cond3 = (all(!is.na(x)))
  cond4 = (all(!is.complex(x)))
  cond5 = (length(dim(x))==2) # only 2-dimensional array is allowed
  cond6 = ((dim(x)[1]!=1)&&(dim(x)[2]!=1))
  
  
  if (cond1&&cond2&&cond3&&cond4&&cond5&&cond6){
    return(TRUE)
  } else {
    stop()
  }
}


n <- 100
p <- 100
mu <- rep(0,p)

# Generate covariance matrix
pre.proportion <- 0.01
pre.value <- 0.3
pre <- matrix(0,p,p)
  ind.d <- sample(length(diag(pre)),round(pre.proportion * length(diag(pre))))
  ind.nd <- sample(sum(lower.tri(pre)),round(pre.proportion * sum(lower.tri(pre))))
  pre[lower.tri(pre)][ind.nd] <- pre.value
  pre <- pre + t(pre)
  diag(pre)[ind.d] <- pre.value
  print(paste("% of nonzero elements: ",sum(pre==0.3)/p^2))
  pre <- pre + (-min(eigen(pre)$value)+0.1^3)*diag(1,p)
cov <- solve(pre)


alps <- seq(0,7,0.01)
gams <- max(2*n, p)^-alps
log.gams <- log(gams / (gams + 1))
FPR.want <- 0.05
n.gen <- 100


given.data <- replicate(2, mvrnorm(n, mu, cov), simplify = "array")
data1 <- given.data[,,1]
data2 <- given.data[,,2]

# Conduct 2-sample mean test with default alpha value
mean2.mxPBF(data1, data2)

# Conduct 2-sample cov test with default alpha value
cov2.mxPBF(data1, data2)

##################################################################################
###################      Empirical FPR method     ################################
##################################################################################

smean <- colMeans(rbind(data1,data2))
svar <- cov(data1)+cov(data2)/2
gen.samples <- replicate(2*n.gen,mvrnorm(n, smean, svar), simplify = "array")


# Find an empirical FPR for mean test
bf.rep <- sapply(1:n.gen, function(l){
  bf <- rep(0,length(alps))
  bf[1] <- mean2.mxPBF(gen.samples[,,l], gen.samples[,,(l+n.gen)], gamma = gams[1])$statistic
  diffs <- cumsum(diff(as.numeric(log.gams)))
  bf[-1] <- sapply(2:length(alps), function(k) bf[1] + diffs[k-1]/2)
  return(bf)
})
empFPR <- apply(bf.rep, 1, function(row) sum(exp(row) > 10)/n.gen)
alps.FPR <- alps[max(which(min(abs(empFPR-FPR.want))==abs(empFPR-FPR.want)))]


# Find an empirical FPR for cov test
bf.rep <- sapply(1:n.gen, function(l){
  bf <- rep(0,length(alps))
  bf[1] <- cov2.mxPBF(gen.samples[,,l], gen.samples[,,(l+n.gen)], gamma = gams[1])$statistic
  diffs <- cumsum(diff(as.numeric(log.gams)))
  bf[-1] <- sapply(2:length(alps), function(k) bf[1] + diffs[k-1]/2)
  return(bf)
})
empFPR <- apply(bf.rep, 1, function(row) sum(exp(row) > 10)/n.gen)
alps.FPR <- alps[max(which(min(abs(empFPR-FPR.want))==abs(empFPR-FPR.want)))]
