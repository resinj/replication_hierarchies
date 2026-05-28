# Load required R packages
library(RcppAlgos)    # comboGeneral(), comboCount(), permuteGeneral()
library(pracma)       # Rank(), nullspace()
library(matlib)       # gaussianElimination
library(RColorBrewer) # brewer.pal()
library(tikzDevice)   # tikz()
library(MASS)         # fractions()


# NOTES:
#   Variable names follow Appendix B of the manuscript
#   k always denotes the number of classes
#   n is the common denominator used in all predicted class probabilities
#     All probabilities are represented as multiples of 1/n!
#   f is a matrix such that f[i,j] = f_j(i) 

# generate matrix with forecast probabilities (multiplied by n)
# on k classes with nonzero class probabilities of the form m/n
# each column represents a forecast
generate_forecasts = function(k,n){
  comb = cbind(comboGeneral(k,n-k,TRUE),matrix(rep(1:k,each = comboCount(k,n-k,TRUE)),ncol = k))
  f = apply(comb,1,table)
  return(f)
}

## Some auxiliary functions
# evaluate the predictive cdfs of forecasts in f at i = 0,1,...,k
# i may be a vector of length N (number of forecasts), 
# in which case the function evaluates the j-th forecast cdf at the j-th value in i
cdf = function(i,f){
  N = ncol(f)
  if(length(i) == 1) i = rep(i,N)
  ifelse(i == rep(0,N), rep(0,N),
         ifelse(i == rep(1,N),f[1,],
                sapply(1:N, function(j) sum(f[1:i[j],j]))))
}

# evaluate the quantile functions of the forecasts at u/n for u between 0 and n
# u may be a vector of length N (number of forecasts), 
# in which case the function evaluates the j-th forecast quantile at the j-th value in u
quant = function(u,f){
  apply(f,2,function(p) 1+sum(cumsum(p) < u))
}

# conditional CDF of PIT evaluated at u/n for each forecast and outcome i
cond_cdf_PIT = function(u,i,f){ 
  ifelse(u < cdf(i-1,f),0,
         ifelse(u >= cdf(i,f),1,
                (u-cdf(i-1,f))/f[i,]))
}

# conditional CDF of the simple PIT (without randomization) evaluated at u/n for each forecast and outcome i
cond_cdf_PIT_simple = function(u,i,f){
  ifelse(u < cdf(i,f),0,1)
}

## linear constraints
# Construct systems of linear equations (sle) corresponding to various notions of calibration
# all functions return a list containing a matrix A, a vector b, 
# and the particular (auto-calibrated) solution f (arranged in a vector) such that Af = b

# linear system for marginal calibration
construct_sle_MC = function(k,n,f){
  N = ncol(f)
  
  # Conditions such that probabilities sum to 1
  A = matrix(rep(diag(N),k),nrow = N)
  b = rep(n,N)
  
  # Conditions for marginal calibration
  A = rbind(A, matrix(rep(diag(k),each = N),nrow = k,byrow = TRUE))
  b = c(b, rowSums(f))
  
  return(list(A = A, b = b, f = as.vector(t(f))))
}

# linear system for probabilistic calibration
construct_sle_PC = function(k,n,f){
  N = ncol(f)
  
  # Conditions such that probabilities sum to 1
  A = matrix(rep(diag(N),k),nrow = N)
  b = rep(n,N)
  
  # Conditions for probabilistic calibration
  matrix_row = function(u,f) as.vector(sapply(1:k,function(i) cond_cdf_PIT(u,i,f)))
  
  u = 1:(n-1)
  A = rbind(A, t(sapply(u,function(u) matrix_row(u,f))))
  b = c(b,u*N)
  
  return(list(A = A, b = b, f = as.vector(t(f))))
}

# linear system for full probabilistic calibration
construct_sle_fPC = function(k,n,f){
  N = ncol(f)
  
  # Conditions such that probabilities sum to 1
  A = matrix(rep(diag(N),k),nrow = N)
  b = rep(n,N)
  
  # Conditions for full probabilistic calibration
  matrix_row = function(u,f,pi) as.vector(sapply(1:k,function(i) cond_cdf_PIT(u,i,f))[,match(1:k,pi)])
  perm = permuteGeneral(k,k)
  
  for(i_pi in 1:nrow(perm)){ 
    u = 1:(n-1)
    pi = perm[i_pi,]
    f_pi = f[pi,]
    A = rbind(A, t(sapply(u,function(u) matrix_row(u,f_pi,pi))))
    b = c(b,u*N)
  }
  
  return(list(A = A, b = b, f = as.vector(t(f))))
}

# linear system for average probabilistic calibration
construct_sle_aPC = function(k,n,f){
  N = ncol(f)
  
  # Conditions such that probabilities sum to 1
  A = matrix(rep(diag(N),k),nrow = N)
  b = rep(n,N)
  
  # Conditions for average probabilistic calibration
  matrix_row = function(u,f,pi) as.vector(sapply(1:k,function(i) cond_cdf_PIT(u,i,f))[,match(1:k,pi)])
  perm = permuteGeneral(k,k)

  A_avg = matrix(0,n-1,k*N)
  b_avg = rep(0,n-1)
  
  for(i_pi in 1:nrow(perm)){ 
    u = 1:(n-1)
    pi = perm[i_pi,]
    f_pi = f[pi,]
    A_avg = A_avg + t(sapply(u,function(u) matrix_row(u,f_pi,pi)))/nrow(perm)
    b_avg = b_avg + u*N/nrow(perm)
  }
  
  return(list(A = rbind(A,A_avg), b = c(b,b_avg), f = as.vector(t(f))))
}

# linear system for conditional (non) exceedance probability calibration
construct_sle_CC = function(k,n,f){
  N = ncol(f)
  
  # Conditions such that probabilities sum to 1
  A = matrix(rep(diag(N),k),nrow = N)
  b = rep(n,N)
  
  # Conditions for conditional (non) exceedance probability calibration
  matrix_row = function(u,f,part){
    vec_row = matrix(0,nrow = N,ncol = k)
    vec_row[part,] = sapply(1:k,function(i) cond_cdf_PIT(u,i,as.matrix(f[,part])))
    as.vector(vec_row)
  }
  
  for(u in seq(1,n-0.5,0.5)){ 
    partition = lapply(1:k,function(i) which(quant(u,f) == i))
    for(i in 1:k){
      part = partition[[i]]
      len = length(part)
      if(len > 0){
        A = rbind(A,matrix_row(u,f,part))
        b = c(b,u*len)
      }
    }
  }
  
  return(list(A = A, b = b, f = as.vector(t(f))))
}

# linear system for full conditional (non) exceedance probability calibration
construct_sle_fCC = function(k,n,f){
  N = ncol(f)
  
  # Conditions such that probabilities sum to 1
  A = matrix(rep(diag(N),k),nrow = N)
  b = rep(n,N)
  
  # Conditions for full conditional (non) exceedance probability calibration
  matrix_row = function(u,f,part,pi){
    vec_row = matrix(0,nrow = N,ncol = k)
    vec_row[part,] = sapply(1:k,function(i) cond_cdf_PIT(u,i,as.matrix(f[,part])))
    as.vector(vec_row[,match(1:k,pi)])
  }
  perm = permuteGeneral(k,k)
  
  for(i_pi in 1:nrow(perm)){
    pi = perm[i_pi,]
    f_pi = f[pi,]
    for(u in seq(1,n-0.5,0.5)){ 
      partition = lapply(1:k,function(i) which(quant(u,f_pi) == i))
      for(i in 1:k){
        part = partition[[i]]
        len = length(part)
        if(len > 0){
          A = rbind(A,matrix_row(u,f_pi,part,pi))
          b = c(b,u*len)
        }
      }
    }
  }
  
  return(list(A = A, b = b, f = as.vector(t(f))))
}

# linear system for threshold calibration
construct_sle_TC = function(k,n,f){
  N = ncol(f)
  
  # Conditions such that probabilities sum to 1
  A = matrix(rep(diag(N),k),nrow = N)
  b = rep(n,N)
  
  # Conditions for threshold calibration
  matrix_row = function(t,f,part){
    vec_row = matrix(0,nrow = N,ncol = k)
    vec_row[part,1:t] = 1
    as.vector(vec_row)
  }
  
  for(t in 1:(k-1)){ 
    partition = lapply(1:(n-1),function(u) which(cdf(t,f) == u))
    for(u in 1:(n-1)){
      part = partition[[u]]
      len = length(part)
      if(len > 0){
        A = rbind(A,matrix_row(t,f,part))
        b = c(b,u*len)
      }
    }
  }
  
  return(list(A = A, b = b, f = as.vector(t(f))))
}

# linear system for full threshold calibration
construct_sle_fTC = function(k,n,f){
  N = ncol(f)
  
  # Conditions such that probabilities sum to 1
  A = matrix(rep(diag(N),k),nrow = N)
  b = rep(n,N)
  
  # Conditions for full threshold calibration
  matrix_row = function(y,f,part,pi){
    vec_row = matrix(0,nrow = N,ncol = k)
    vec_row[part,1:y] = 1
    as.vector(vec_row[,match(1:k,pi)])
  }
  
  perm = permuteGeneral(k,k)
  
  for(i_pi in 1:nrow(perm)){
    pi = perm[i_pi,]
    f_pi = f[pi,]
    for(t in 1:(k-1)){ 
      partition = lapply(1:(n-1),function(u) which(cdf(t,f_pi) == u))
      for(u in 1:(n-1)){
        part = partition[[u]]
        len = length(part)
        if(len > 0){
          A = rbind(A,matrix_row(t,f_pi,part,pi))
          b = c(b,u*len)
        }
      }
    }
  }
  
  return(list(A = A, b = b, f = as.vector(t(f))))
}

# linear system for double probabilistic calibration
construct_sle_DC = function(k,n,f){
  N = ncol(f)
  
  # Conditions such that probabilities sum to 1
  A = matrix(rep(diag(N),k),nrow = N)
  b = rep(n,N)
  
  # Conditions for double probabilistic calibration
  matrix_row = function(u,f) as.vector(sapply(1:k,function(i) cond_cdf_PIT_simple(u,i,f)))
  
  u = 1:(n-1)
  A = rbind(A,t(sapply(u,function(u) matrix_row(u,f))))
  b = c(b,sapply(u,function(u) sum(sapply(1:k,function(i) f[i,]*cond_cdf_PIT_simple(u,i,f)))))
  
  return(list(A = A, b = b, f = as.vector(t(f))))
}

# linear system for class-wise calibration
construct_sle_CwC = function(k,n,f){
  N = ncol(f)
  
  # Conditions such that probabilities sum to 1
  A = matrix(rep(diag(N),k),nrow = N)
  b = rep(n,N)
  
  # Conditions for class-wise calibration
  matrix_row = function(i,f,part){
    vec_row = matrix(0,nrow = N,ncol = k)
    vec_row[part,i] = 1
    as.vector(vec_row)
  }
  
  for(i in 1:k){ 
    partition = lapply(1:(n+1-k),function(u) which(f[i,] == u))
    for(u in 1:(n+1-k)){
      part = partition[[u]]
      len = length(part)
      if(len > 0){
        A = rbind(A,matrix_row(i,f,part))
        b = c(b,u*len)
      }
    }
  }
  
  return(list(A = A, b = b, f = as.vector(t(f))))
}

# function that tests whether the vector g (or f if a solution g is not provided)
# satisfies Ag = b
test_sle = function(sle, g = NULL){
  if(is.null(g)) g = sle$f
  else g = as.vector(t(g))
  
  all(abs(sle$A %*% g - sle$b) < 10^-8) # avoid numerical instabilities
}

# function that tests which of the notions is satisfied by the specified example
test_all = function(k,n,f,g){
  sle_CwC = construct_sle_CwC(k,n,f)
  CwC = test_sle(sle_CwC,g)
  
  sle_MC = construct_sle_MC(k,n,f)
  MC = test_sle(sle_MC,g)
  
  sle_fCC = construct_sle_fCC(k,n,f)
  fCC = test_sle(sle_fCC,g)
  
  # sle_CC = construct_sle_CC(k,n,f)
  # CC = test_sle(sle_CC,g)
  
  sle_fPC = construct_sle_fPC(k,n,f)
  fPC = test_sle(sle_fPC,g)
  
  # sle_PC = construct_sle_PC(k,n,f)
  # PC = test_sle(sle_PC,g)
  
  sle_aPC = construct_sle_aPC(k,n,f)
  aPC = test_sle(sle_aPC,g)
  
  sle_fTC = construct_sle_fTC(k,n,f)
  fTC = test_sle(sle_fTC,g)
  
  # sle_TC = construct_sle_TC(k,n,f)
  # TC = test_sle(sle_TC,g)
  
  # partial versions
  perm = permuteGeneral(k,k)

  pCC = FALSE
  pPC = FALSE
  pTC = FALSE
  for(i_pi in 1:nrow(perm)){
    pi = perm[i_pi,]
    f_pi = f[pi,]
    g_pi = g[pi,]
    
    sle_CC = construct_sle_CC(k,n,f_pi)
    pCC = pCC || test_sle(sle_CC,g_pi)
    
    sle_PC = construct_sle_PC(k,n,f_pi)
    pPC = pPC || test_sle(sle_PC,g_pi)
    
    sle_TC = construct_sle_TC(k,n,f_pi)
    pTC = pTC || test_sle(sle_TC,g_pi)
  }

  cbind(CwC,MC,fCC,pCC,fPC,pPC,aPC,fTC,pTC)
}

# combine two linear systems
combine_sle = function(sle1,sle2){
  if(!all(sle1$f == sle1$f)) stop("incompatible systems")
  sle = list(A = rbind(sle1$A,sle2$A), b = c(sle1$b,sle2$b), f = sle1$f)
  return(sle)
}

# find the homogenous solution for a linear system
solve_hom = function(sle){
  # print(paste0("Rank = ",Rank(sle$A)))
  hom_sol = nullspace(sle$A)
  # find a nice basis for the kernel
  if(dim(hom_sol)[[2]] > 1) hom_sol = t(gaussianElimination(t(hom_sol)))
  else hom_sol = hom_sol/hom_sol[head(which(abs(hom_sol) > 10^-12),1)]
  return(hom_sol)
}

## Visualizations
# unconditional coverage plot 
plot_UQC = function(k,n,f,g){
  N = ncol(f)
  
  coverage_lower = function(u){
    sapply(u, function(u) cdf(quant(u,f) - 1,g))
  }
  coverage_upper = function(u){
    sapply(u, function(u) cdf(quant(u,f),g))
  }
  par(mfrow = c(1,1), lwd = 2,mar = c(2.5,2.5,0.1,0.1),mgp = c(1.4,0.4,0))
  plot(NULL,xlim = c(0,1),ylim = c(0,1),xlab = expression(alpha),ylab = "coverage")
  cov_l = coverage_lower(0:n)/n
  cov_u = coverage_upper(0:n)/n
  abline(0,1,col = "gray")
  # for(j in 1:N){
  #   lines(0:n/n, cov_l[j,],type = "S",col = colors[j],lty = "dashed")
  #   lines(0:n/n, cov_u[j,],type = "S",col = colors[j],lty = "longdash")
  # }
  lines(0:n/n, colMeans(cov_l),type = "S",lty = "dashed")
  lines(0:n/n, colMeans(cov_u),type = "S",lty = "longdash")
  
  legend("topleft",legend = c("upper","lower"),lty = c("longdash","dashed"))
}

# conditional coverage and quantile functions for each forecast
plot_QC = function(k,n,f,g){
  N = ncol(f)
  
  colors = brewer.pal(N,"Set2")
  colors_adj = adjustcolor(colors, alpha.f = 0.6)

  coverage_lower = function(u){
    sapply(u, function(u) cdf(quant(u,f) - 1,g))
  }
  coverage_upper = function(u){
    sapply(u, function(u) cdf(quant(u,f),g))
  }
  # Coverages
  par(mfrow = c(1,1), lwd = 2,mar = c(2.5,2.5,0.1,0.1),mgp = c(1.4,0.4,0))
  plot(NULL,xlim = c(0,1),ylim = c(0,1),xlab = expression(alpha),ylab = "coverage")
  cov_l = coverage_lower(0:n)/n
  cov_u = coverage_upper(0:n)/n
  abline(0,1,col = "gray")
  for(j in 1:N){
    lines(0:n/n, cov_l[j,],type = "S",col = colors_adj[j],lty = "dashed")
    lines(0:n/n, cov_u[j,],type = "S",col = colors_adj[j],lty = "longdash")
  }
  lines(0:n/n, colMeans(cov_l),type = "S",lty = "dashed")
  lines(0:n/n, colMeans(cov_u),type = "S",lty = "longdash")
  
  legend("topleft",legend = c("upper","lower"),lty = c("longdash","dashed"))
  legend("bottomright",legend = as.expression(sapply(1:3, \(j) bquote(F[.(j)]))),col = colors_adj,lty = 1)
  # Quantiles
  plot(NULL,xlim = c(0,1),ylim = c(1,k),xlab = expression(alpha),ylab = "quantile")
  quantiles = sapply(0:n, function(u) quant(u,f))
  for(j in 1:N){
    lines(0:n/n, quantiles[j,],type = "S",col = colors[j],lty = j)
  }
  legend("bottomright",legend = as.expression(sapply(1:3, \(j) bquote(F[.(j)]))),col = colors_adj[1:N],lty = 1:N)
}
