#' @name comire.gibbs
#' 
#' @export comire.gibbs
#'
#' @title Gibbs sampler for CoMiRe model
#' 
#' @usage comire.gibbs(y, x, z = NULL, family = 'continuous', 
#'        grid = NULL, mcmc, prior, 
#'        state = NULL, seed, max.x = max(x), z.val = NULL)
#' 
#' @description Posterior inference via Gibbs sampler for CoMiRe model
#' 
#' @details 
#' The function fit a convex mixture regression (\code{CoMiRe}) model (Canale, Durante, Dunson, 2018) via Gibbs sampler. 
#' For continuous outcome \eqn{y \in \mathcal{Y}}{y}, adverse esposure level \eqn{x \in \mathcal{X}}{x} and no confunding 
#' variables, one can set \code{family = 'continuous'} and \code{z = NULL} and fit model
#' \cr
#' \eqn{ f_x(y) = \{1-\beta(x)\} \sum_{h=1}^{H}\nu_{0h} \phi(y; \theta_{0h}, \tau_{0h}^{-1}) + \beta(x) \phi(y; \theta_{\infty}, \tau_{\infty}^{-1})}{ f(y | x) = {1-\beta(x)} \sum \nu_0h \phi(y; \theta_0h, 1/(\tau_0h)) + \beta(x) \phi(y; \theta_{inf}, 1/(\tau_{inf}))} ;\cr
#' \cr
#' where \eqn{\beta(x) = \sum_{j=1}^{J} \omega_j \psi_j(x), x\ge0,}{\beta(x) = \sum \omega_j \psi_j(x), x\ge0,} 
#' is a a monotone nondecreasing interpolation function, constrained between 0 and 1 and  \eqn{\psi_1,...,\psi_J} are monotone nondecreasing I-splines basis. 
#' \cr
#' If \eqn{p \ge 1} confounding covariates \eqn{z \in \mathcal{Z}}{z} are available, passing the argument \code{z} 
#' the function fits model\cr
#' \cr
#' \eqn{f(y; x,z) = \{1-\beta(x)\} f_0(y;z) + \beta(x) f_\infty(y;z)}{f(y | x,z) = {1-\beta(x)} f_0(y | z) + \beta(x) f_{inf}(y | z)} ;\cr
#' \cr
#' where: \cr
#' \eqn{f_0(y;z)= \sum_{h=1}^{H} \nu_{0h} \phi(y;\theta_{0h}+z^\mathsf{T}\gamma,\tau_{0h}^{-1})}{f_0(y | z)= \sum \nu_0h \phi(y; \theta_0h + z' \gamma, 1/(\tau_0h) )}, and 
#' \eqn{f_\infty(y;z)= \phi(y;\theta_\infty+ z^\mathsf{T}\gamma,\tau_{\infty}^{-1})}{f_{inf}(y | z)= \phi(y; \theta_{inf} + z' \gamma, 1/(\tau_{inf}))}. \cr
#' \cr
#' Finally, if \eqn{y} is a binary response, one can set \code{family = 'binary'} and fit model \cr
#' \cr
#' \eqn{p_x(y) = (\pi_x)^y (1 - \pi_x)^{1-y}}{p(y | x) = (\pi_x)^y (1 - \pi_x)^(1-y)} ; \cr
#' \cr
#' where \eqn{\pi_x = P(Y=1 | x)} is 
#' \eqn{\pi_x = \{1-\beta(x)\} \pi_0 + \beta(x) \pi_\infty}{\pi_x = {1-\beta(x)} \pi_0 + \beta(x) \pi_inf}.
#' 
#'
#' @param y numeric vector for the response:
#' when \code{family="continuous"} \code{y} must be a numeric vector; if \code{family="binary"} \code{y} must assume values \code{0} or \code{1}.
#' 
#' @param x numeric vector for the covariate relative to the dose of exposure.
#' 
#' @param z numeric vector for the confunders; a vector if there is only one 
#' confounder or a matrix for two or more confunders
#' 
#' @param family type of \code{y}. This can be \code{"continuous"} or \code{"binary"}. Default  \code{"continuous"}.
#' 
#' @param grid a list giving the parameters for plotting the posterior mean density and the posterior mean \eqn{\beta(x)} over finite grids
#'  if \code{family="continuous"} and \code{z=NULL}. It must include the following values:
#'  \itemize{
#'  \item \code{grids}, logical value (if \code{TRUE} the provided grids are used, otherwise standard grids are automatically used); 
#'  \item \code{xgrid} and \code{ygrid}, numerical vectors with the actual values of the grid for y and x.
#'  }
#'
#' @param mcmc a list giving the MCMC parameters. It must include the following integers: \code{nb} giving the number of burn-in iterations, \code{nrep} giving the total number of iterations, \code{thin} giving the thinning interval, \code{ndisplay} giving the multiple of iterations to be displayed on screen while the algorithm is running (a message will be printed every \code{ndisplay} iterations).
#'
#' @param prior a list containing the values of the hyperparameters. 
#' 
#' If \code{family = "continuous"}, it must include the following values: 
#' \itemize{
#' \item \code{mu.theta}, the prior mean \eqn{\mu_\theta} for each location parameter \eqn{\theta_{0h}}{\theta_0h} and \eqn{\theta_1}, 
#' \item \code{k.theta}, the prior variance \eqn{k_\theta} for each location paramter \eqn{\theta_{0h}}{\theta_0h} and \eqn{\theta_1}, 
#' \item \code{mu.gamma} (if \code{p} confounding covariates are included in the model) a \code{p}-dimentional vector of prior means \eqn{\mu_\gamma}{\mu_gamma} of the parameters \eqn{\gamma} corresponding to the confounders,
#' \item \code{k.gamma}, the prior variance \eqn{k_\gamma}{k_gamma} for parameter corresponding to the confounding covariate (if \code{p=1}) or \code{sigma.gamma} (if \code{p>1}), that is the covariance matrix \eqn{\Sigma_\gamma}{\Sigma_gamma} for the parameters corresponding to the \code{p} confounding covariates; this must be a symmetric positive definite matrix.
#' \item \code{eta}, numeric vector of size \code{J} for the Dirichlet prior on the beta basis weights, 
#' \item \code{alpha}, prior for the mixture weights,
#' \item \code{a} and \code{b}, prior scale and shape parameter for the gamma distribution of each precision parameter, 
#' \item \code{J}, parameter controlling the number of elements of the I-spline basis,
#' \item \code{H}, total number of components in the mixture at \eqn{x_0}.
#' }
#' 
#' If \code{family="binary"} it must include the following values: 
#' \itemize{
#' \item \code{eta}, numeric vector of size \code{J} for the Dirichlet prior on the beta basis weights, 
#' \item \code{a.pi0} and \code{b.pi0}, the prior parameters of the prior beta distribution for \eqn{\pi_0},
#' \item \code{J}, parameter controlling the number of elements of the Ispline basis.
#' }
#' 
#' @param state if \code{family="continuous"}, a list giving the current value of the parameters. This list is used if the current analysis is the continuation of a previous analysis or if we want to start the MCMC algorithm from some particular value of the parameters.
#' 
#' @param seed seed for random initialization.
#' 
#' @param max.x maximum value allowed for \code{x}.
#' 
#' @param z.val optional numeric vector containing a fixed value of interest for each of the confounding covariates to be used for the plots. Default value is \code{mean(z)} for numeric covariates or the mode for factorial covariates.
#'
#' @return An object of the class \code{classCoMiRe}, i.e. a list of arguments for generating posterior output. It contains:
#' \itemize{
#' \item{\code{call}}{the model formula}
#' \item{\code{post.means}}{ a list containing the posterior mean density beta over the grid, of all the mixture parameters and, 
#' if \code{family = "continuous"} and \code{z = NULL}, of \eqn{f_0} and \eqn{f_{inf}} over the \code{y.grid}.}
#' \item{\code{ci}}{ a list containing the 95\% credible intervals for all the quantities stored in \code{post.means}.}
#' \item{\code{mcmc}}{ a list containing all the MCMC chains.}
#' \item{\code{z}}{ the same of the input}
#' \item{\code{z.val}}{ the same of the input}
#' \item{\code{f0,f1}}{ MCMC replicates of the density in the two extremes (only if \code{family = 'continuous'})}
#' \item{\code{nrep,nb}}{ the same values of the list \code{mcmc} in the input arguments}
#' \item{\code{bin}}{ logical, equal to \code{TRUE} if \code{family = 'binary'}}
#' \item{\code{univariate}}{ logical, equal to \code{TRUE} if \code{z} is null or a vector}
#' } 
#' 
#' @references Canale, A., Durante, D., and Dunson, D. (2018), Convex Mixture Regression for Quantitative Risk Assessment, Biometrics, 74, 1331-1340
#' @references Galtarossa, L., Canale, A., (2019), A Convex Mixture Model for Binomial Regression, Book of Short Papers SIS 2019
#'
#' @author Antonio Canale [aut, cre], Daniele Durante [ctb], Arianna Falcioni [aut], Luisa Galtarossa [aut], Tommaso Rigon [ctb]
#'
#' @examples{
#' \dontrun{
#' 
#' attach(CPP)
#' 
#' n <- NROW(CPP)
#' J <- H <- 10
#' 
#' mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#' 
#' 
#' ## 1. continuous case ##
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, eta=rep(1, J)/J, 
#'               alpha=rep(1,H)/H, a=2, b=2, J=J, H=H)
#' fit0 <- comire.gibbs(gestage, dde, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=5, max.x=180)
#' 
#' 
#' ## 1.1 only one confounding covariate  ##
#' 
#' mage_std <- scale(mage, center = T, scale = T) 
#' 
#' prior <- list(mu.theta=mean(gestage), k.theta=10, mu.gamma=0, k.gamma=10, 
#'               eta=rep(1, J)/J, alpha=1/H, a=2, b=2, H=H, J=J)
#' fit1 <- comire.gibbs(gestage, dde, mage_std, family="continuous", 
#'               mcmc=mcmc, prior=prior, seed=5, max.x=180)
#' 
#' 
#' ## 1.2 more than one confounding covariates  ##
#' 
#' z <- cbind(mage, mbmi, mbmi^2, sei)
#' z <- scale(z, center = T, scale = T) 
#' z <- as.matrix(cbind(z, CPP$smoke))
#' mod<- lm(gestage ~ dde + )
#' 
#' prior <- list(mu.theta=mod$coefficients[1], k.theta=10, 
#'               mu.gamma=mod$coefficients[-c(1,2)], sigma.gamma=diag(rep(10,5)), 
#'               eta=rep(1, J)/J, alpha=1/H, a=2, b=2, H=H, J=J)
#' fit2 <- comire.gibbs(y, x, z, family="continuous", 
#'                      mcmc=mcmc, prior=prior, seed=5)
#' 
#' 
#' ## 2. binary case ##
#' 
#' premature <- as.numeric(y<37)
#' 
#' prior <- list(pi0=mean(gestage), eta=rep(1, J)/J, 
#'              a.pi0=27, b.pi0=360, J=J)
#' fit_binary<- comire.gibbs(premature, dde, family="binary", 
#'                           mcmc=mcmc, prior=prior, seed=5, max.x=180)
#' 
#' }
#' }
#' 
#' 
#' 

comire.gibbs<- function(y, x, z = NULL, family = 'continuous', 
                        grid = NULL, mcmc, prior, state=NULL, 
                        seed, max.x = max(x), z.val=NULL){
    
    # check family
    if(is.null(family)) stop ("Missing family parameter")
  
    # family="continuous"|"binary"
    if(!(family=="continuous"|family=="binary")) stop("Wrong family specification")
    
    # controllo valori 0|1
    if(family=="binary" & any((levels(factor(y))==c(0,1))==F)) stop("Values of y must be 0 or 1")
    
    # per caso binarye non posso usare covariate
    if(family=="binary" & !is.null(z)) stop("Confounders not (yet) implemented if family = 'binary'")

    if(family=="binary" & is.null(z)) {
      cat("CoMiRe model fit via Gibbs Sampler\n")
      cat("Family: binary\n")
      call <- paste(c(deparse(substitute(y)), " ~ . | beta(", deparse(substitute(x)), ")"), collapse="")  
      res <- comire.gibbs.binary(y, x, mcmc, prior, seed, max.x=max(x))
      as.classCoMiRe(call= call, out = res, nrep= mcmc$nrep, nb= mcmc$nb, bin = TRUE)
    }

    else if(family=="continuous" & is.null(z)) {
      cat("CoMiRe model fit via Gibbs Sampler\n")
      cat("Family: continuous\n")
      call <- paste(c(deparse(substitute(y)), " ~ . | beta(", 
                      deparse(substitute(x)), ")"), collapse="")  
      res <- comire.gibbs.continuous(y, x, grid=grid, mcmc, prior, state=state, 
                                     seed, max.x=max(x))
      as.classCoMiRe(call= call, out = res, nrep= mcmc$nrep, nb= mcmc$nb)
    }
  
    else if(family=="continuous" & !is.null(z) & (is.null(ncol(z)) | ncol(z)==1) ) {
      cat("CoMiRe model fit via Gibbs Sampler\n")
      cat("Family: continuous \n")
      call <- paste(c(deparse(substitute(y)), " ~ ", 
                      deparse(substitute(z)), " | beta(", 
                      deparse(substitute(x)), ")"), collapse="")  
      res <- comire.gibbs.continuous.confunder(y, x, z, grid=grid, mcmc, prior, 
                                               state=state, seed, max.x=max(x),
                                               z.val=z.val)
      as.classCoMiRe(call= call, out = res$out,  z = z, z.val = res$z.val, 
                     nrep = mcmc$nrep, nb = mcmc$nb)
    }
  
    else if(family=="continuous" & !is.null(z) & !is.null(ncol(z))) {
      cat("CoMiRe model fit via Gibbs Sampler\n")
      cat("Family: continuous \n")
      call <- paste(c(deparse(substitute(y)), " ~ ",
                      paste(attr(z, "dimnames")[[2]], collapse =" + "), 
                      " | beta(", deparse(substitute(x)), ")"), collapse="")
      res <- comire.gibbs.continuous.confunders(y, x, z, grid=grid, mcmc, prior, 
                                                state=state, seed, max.x=max(x), 
                                                z.val=z.val)
      as.classCoMiRe(call = call, out = res$out, z = z, z.val = res$z.val, 
                     nrep = mcmc$nrep, nb = mcmc$nb, 
                     univariate = FALSE)
      
    }

}
