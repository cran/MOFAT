#' @title
#' MOFAT
#'
#' @description
#' This function can be used for generating MOFAT designs.
#' 
#' @details 
#' The \code{mofat} function generates the MOFAT design 
#' for a given number of factors (\eqn{p\ge 2}) and 
#' number of base runs (\eqn{l \ge 3}). The total number of runs in the MOFAT design will be \eqn{l(p+1)}. 
#' A MOFAT design can be viewed as an optimized version of Morris screening design (Morris 1991) by exploiting 
#' its connections with the Monte Carlo-based design of Sobol' (1993).
#' Please see Xiao et al. (2022) for details.
#' 
#' Three choices for the \code{method} are given: "uniform", "projection", and "best". Option "uniform" gives \code{l} equally-spaced levels 
#' for the entire design, which are also balanced. "projection" option adjusts the levels of the two base matrices A and B such that 
#' there are \eqn{2l} or \eqn{2l-1} levels in the design depending on \code{l} is even or odd. Option "best" (default) chooses the best 
#' among the first two options using maximin distance criterion.
#' 
#' @author Qian Xiao and V. Roshan Joseph
#' 
#' @importFrom SLHD maximinSLHD
#' @importFrom stats dist
#' 
#' @param p number of factors
#' @param l number of base runs 
#' @param method choose among "uniform", "projection", and "best" 
#' 
#' @return 
#' \item{design}{MOFAT design}
#' 
#' @references 
#' 
#' Morris, M. D. (1991), “Factorial sampling plans for preliminary computational experiments,”
#' Technometrics, 33, 161–174.
#'
#' Sobol’, I. M. (1993), “On sensitivity estimation for nonlinear mathematical models,” Mathematical
#' Modeling and Computational Experiments, 1, 407–414.
#'
#' Xiao, Q., Joseph, V. R., and Ray, D. M. (2022). “Maximum One-Factor-At-A-Time  Designs for Screening in Computer Experiments”.   
#' Technometrics, to appear.
#'
#' @export
#' 
#' @examples
#' #MOFAT with three base runs
#' mofat(p=10, l=3, method="uniform")
#' mofat(p=10, l=3, method="projection")
#' 
#' #MOFAT with five base runs
#' mofat(p=10,l=5)
#' dim(mofat(p=125,l=5))

#library(SLHD)

###########################################
####MOFAT design with l levels

mofat <- function(p, l, method="best") #p is the factor size; l is the level size.
{
  #number of runs needed
  n=l*(p+1)
  
  ##################function needed
  #function for level permutation
  permut <- function(n){
    if(n==1){
      return(matrix(1))
    } else {
      sp <- permut(n-1)
      p <- nrow(sp)
      A <- matrix(nrow=n*p,ncol=n)
      for(i in 1:n){
        A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
      }
      return(A)
    }
  }
  
  #generate B matrix based on A matrix
  B.gen = function(A, l)
  {
    B = A
    for (j in 1:dim(A)[2]) {
      B[,j] =  (rank(A[,j] -.5+as.numeric(A[,j]<.5))-0.5)/l
    }
    return(B)
  }
  
  #generate design based on A and B
  sobol_design = function(A,B)
  {
    d = dim(A)[2]
    res = A
    for (i in 1:d) {
      interm = A
      interm[,i] = B[,i]
      res = rbind(res, interm)
    }
    return(res)
  }
  
  adjust=function(A,B){
    A0=A
    B0=B
    A[A<.5]=A[A<.5]-.25/l
    A[A>.5]=A[A>.5]+.25/l
    B[B<.5]=B[B<.5]+.25/l
    B[B>.5]=B[B>.5]-.25/l
    D1=sobol_design(A,B)
    d1=min(dist(D1,method = "manhattan"))
    
    A=A0
    B=B0
    A[A<.5]=A[A<.5]+.25/l
    A[A>.5]=A[A>.5]-.25/l
    B[B<.5]=B[B<.5]-.25/l
    B[B>.5]=B[B>.5]+.25/l
    D2=sobol_design(A,B)
    d2=min(dist(D2,method = "manhattan"))
    if(d2<d1) design=D1 else design=D2
    return(design)
  }
  
  choose_design=function(A,B){
    design = sobol_design(A,B)
    d0=min(dist(design,method = "manhattan"))
    if(method=="projection") design=adjust(A,B)
    if(method=="best") {
      design_adj=adjust(A,B)
      d_adj=min(dist(design_adj,method = "manhattan"))
      if(d_adj>d0) design=design_adj
    }
    return(design)
  }
  ###Mofat with l=3 case:
  if (l == 3)
  {
    ###when the number of factors k is no larger than l!
    if (p <= 6)
    {
      #creation of matrix P_l
      allperm = t(permut(l))
      A = allperm[, 1:p]
      ###standardize to (0,1) range.
      A = (A-0.5)/l
      B = B.gen(A,l)
      design=choose_design(A,B)
      return(design)
    } else {
      #creation of matrix P_l
      allperm = t(permut(l))
      #find the best A from Algebraic construction
      if(p%%6 == 0)
      {
        A = matrix(rep(allperm, floor(p/6)), l, floor(p/6)*dim(allperm)[2])
      } else {
        A = matrix(rep(allperm, floor(p/6)), l, floor(p/6)*dim(allperm)[2])
        A = cbind(A, allperm[,1:(p%%6)])
      }
      ###standardize to (0,1) range.
      A = (A-0.5)/l
      B = B.gen(A,l)
      design=choose_design(A,B)
      return(design)
    } 
  } else {
    ###when the number of factors k is no larger than l!
    if (p <= factorial(l))
    {
      A = maximinSLHD(t = 1, m = l, k = p, nstarts = 100)$StandDesign
      #A=LHD(l,p)
      
      #matrix B
      B = B.gen(A,l)
      design=choose_design(A,B)
      return(design)
      
    } else {
      #creation of matrix P_l
      allperm = t(permut(l))
      #column number in allperm
      cn = dim(allperm)[2]
      #rescale allperm
      allperm = (allperm-0.5)/l
      
      #A' matrix
      A_p = maximinSLHD(t = 1, m = l, k = (p%%cn), nstarts = 100)$StandDesign
      #A matrix
      A = matrix(rep(allperm, floor(p/cn)), l, floor(p/cn)*cn)
      A = cbind(A, A_p)
      #matrix B
      B = B.gen(A,l)
      design=choose_design(A,B)
      return(design)
    }
  }
}

#################################
#' @title
#' Screening measures
#'
#' @description
#' This function can be used for computing screening measures.
#' 
#' @details 
#' The \code{measure} function computes the screening measures such as the total Sobol' indices (Sobol' 1993) 
#' and \eqn{\mu^*} measure of Campolongo et al. (2007). The design matrix should have the Sobol' design structure. 
#' Please see Xiao et al. (2022) for details.
#' 
#' @author Qian Xiao and V. Roshan Joseph
#' 
#' @importFrom stats var
#' 
#' @param design design matrix, which should have the Sobol' design structure
#' @param y response vector
#' 
#' @return 
#' \item{t}{Total Sobol' index}
#' \item{mustar}{\eqn{\mu^*} measure}
#' 
#' @references 
#' 
#'Sobol’, I. M. (1993), “On sensitivity estimation for nonlinear mathematical models,” Mathematical
#'Modeling and Computational Experiments, 1, 407–414.
#'
#'Campolongo, F., Cariboni, J., and Saltelli, A. (2007), “An effective screening design for
#'sensitivity analysis of large models,” Environmental modelling and software, 22, 1509–1518.
#'
#'Xiao, Q., Joseph, V. R., and Ray, D. M. (2022).   
#' “Maximum One-Factor-At-A-Time  Designs for Screening in Computer Experiments”. Technometrics, to appear.
#'
#' @export
#' 
#' @examples
#' #Friedman function
#' fun <- function (X)
#' {
#'  Y <- 10*sin(pi*X[1]*X[2]) + 20*(X[3] - 0.5)^2 + 10*X[4] + 5*X[5]
#'  return(Y)
#' }
#' design = mofat(p=10, l=3)
#' y = apply(design, 1, fun)
#' 
#' #Screening measures
#' measure(design, y)


measure <- function(design, y) 
{
  design = as.matrix(design)
  y = as.vector(y)
  #number of data
  n = dim(design)[1]
  #number of factor
  p = dim(design)[2]
  #number of levels
  l = n/(p+1)
  
  #calculate V_tot
  V_tot = numeric(p)
  #response for matrix A
  y0 = y[1:l]
  for (i in 1:p) {
    V_tot[i] = sum((y0 - y[(i*l+1):(i*l+l)])^2) / (2*l)
  }
  
  #calculate Total Sobol' index
  TS = V_tot/var(y)
  
  #calculate mu star value
  mu_star = numeric(p)
  #matrix A
  A = design[1:l, ]
  #response for matrix A
  y0 = y[1:l]
  for (i in 1:p) {
    #Corresponding matrix C
    C = design[(i*l+1):(i*l+l), ]
    mu_star[i] = sum(abs((y0 - y[(i*l+1):(i*l+l)])/(A[,i] - C[,i])))/l
  }
  return(list(t = TS, mustar = mu_star))
}
