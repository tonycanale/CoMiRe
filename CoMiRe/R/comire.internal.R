#' Internal functions of the package
#' @name comire.internal
#' 
#' 
#' @keywords internal
#'


#-----------------
#Internal functions comire.continuous

#' pssq_gaussian 
#' 
#' @param index index
#' @param data data
#' @param cluster cluster
#' @param locations locations
#' 
pssq_gaussian <- function(index, data, cluster, locations) sum((data[cluster==index] - locations[index])^2)

#-----------------
#Internal functions comire continuous with confounders

#' labelling_b_uni 
#' 
#' @param i i
#' @param w w
#' @param phi phi
#' @param f0i f0i
#' @param f1i f1i
#' 
labelling_b_uni <- function(i, w, phi, f0i, f1i)
{
  probs <- w*((1-phi[i,])*f0i[i] + phi[i,]*f1i[i])
  if(any(probs<0)) probs[probs<0]=0
  sample(1:ncol(phi), 1, prob=probs)
}

#

#' labelling_c_uni 
#' 
#' @param i i
#' @param y y
#' @param z z
#' @param nu nu
#' @param theta theta
#' @param tau tau
#' @param ga ga
#' 
labelling_c_uni <- function(i, y, z, nu, theta, tau, ga)
{
  probs <- nu*stats::dnorm(y[i], theta+z[i]*ga , sqrt(1/tau))
  if(any(probs<0)) probs[probs<0]=0
  sample(1:length(nu), 1, prob=probs)
}

#

#' mixdensity_uni 
#' 
#' @param i i
#' @param y y
#' @param z z
#' @param nu nu
#' @param theta theta
#' @param tau tau
#' @param ga ga
#' 
mixdensity_uni <- function(i, y, z, nu, theta, tau, ga)
{
  kernels <- stats::dnorm(y[i], (theta+z[i]*ga) , sqrt(1/tau))
  fji <- sum(nu * kernels)
  fji
}

#

#' pssq_uni 
#' 
#' @param index index
#' @param y y
#' @param z z
#' @param cluster cluster
#' @param theta theta
#' @param gamma gamma
#' 
pssq_uni <- function(index, y, z, cluster, theta, gamma) {
  sum((y[cluster==index] - 
         (theta[index]+z[cluster==index]*gamma))^2)}

#

#' psdp_uni (posterior double product)
#' 
#' @param index index
#' @param y y
#' @param z z
#' @param cluster cluster
#'
psdp_uni <- function(index, y, z, cluster){
  sum(y[cluster==index]*z[cluster==index])
}

#

#' labelling_b_multi
#' 
#' @param i i
#' @param w w
#' @param phi phi
#' @param f0i f0i
#' @param f1i f1i
#'
labelling_b_multi <- function(i, w, phi, f0i, f1i)
{
  probs <- w*((1-phi[i,])*f0i[i] + phi[i,]*f1i[i])
  if(any(probs<0)) probs[probs<0]=0
  sample(1:ncol(phi), 1, prob=probs)
}

#

#' labelling_c_multi
#' 
#' @param i i
#' @param y y
#' @param z z
#' @param nu nu
#' @param theta theta
#' @param tau tau
#' @param ga ga
#'
labelling_c_multi <- function(i, y, z, nu, theta, tau, ga)
{
  probs <- nu*stats::dnorm(y[i], theta+as.vector(crossprod(z[i,],ga)) , sqrt(1/tau))
  if(any(probs<0)) probs[probs<0]=0
  sample(1:length(nu), 1, prob=probs)
}

#

#' mixdensity_multi
#' 
#' @param i i
#' @param y y
#' @param z z
#' @param nu nu
#' @param theta theta
#' @param tau tau
#' @param ga ga
#'
mixdensity_multi <- function(i, y, z, nu, theta, tau, ga) 
{
  kernels <- stats::dnorm(y[i], theta+as.vector(crossprod(z[i,],ga)) , sqrt(1/tau))
  fji <- sum(nu * kernels)
  fji
}

#

#' pssq_multi
#' 
#' @param y y
#' @param z z
#' @param times times
#' @param gamma gamma
#' @param theta theta
#'
pssq_multi <- function(y,times,z,gamma,theta){
  H <- length(theta)
  sapply(1:H, function(x) 
    crossprod( y-(theta[x]*rep(1,times)+crossprod(t(z),gamma)) ))
}

#

#' psdp_multi
#' 
#' @param index index
#' @param y y
#' @param z z
#' @param cluster cluster
#'
psdp_multi <- function(index, y, z, cluster){
  sum(y[cluster==index]*z[cluster==index])
}

