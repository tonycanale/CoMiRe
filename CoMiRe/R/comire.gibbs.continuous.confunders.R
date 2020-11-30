comire.gibbs.continuous.confunders <- function(y, x, z, grid=NULL, mcmc, prior, 
                            state=NULL, seed=5, max.x=ceiling(max(x)),
                            z.val=NULL, verbose = TRUE){
  
  
  if(is.null(z.val)){
    z.val <- apply(z, 2, function(x) if(length(table(x))>2){mean(x)} 
                   else{unique(x)[which.max(tabulate(match(x,unique(x))))]})
  }
  
  # internal working variables
  n <- length(y)
  p <- ncol(z)
  H <- prior$H 
  J <- prior$J
  print_now <- c(mcmc$nb + 1:mcmc$ndisplay*(mcmc$nrep)/mcmc$ndisplay)
  
  # create the objects to store the MCMC samples
  w <- matrix(NA, mcmc$nrep+mcmc$nb, J)
  nu0 <- matrix(NA, mcmc$nrep+mcmc$nb, H)         
  th0 <- tau0 <- matrix(NA, mcmc$nrep+mcmc$nb, H)
  nu1 <- th1 <- tau1 <- rep(NA, mcmc$nrep+mcmc$nb)
  ga <- matrix(NA, mcmc$nrep+mcmc$nb, p)
  
  # initialize each quantity
  
  ## paramters of the model
  if(is.null(state))
  {
    set.seed(seed)
    w[1,] <- prior$eta/sum(prior$eta)
    nu0[1,] <- rep(1/H, H)
    tau0[1,] <- rgamma(H, prior$a, prior$b)
    th0[1,] <- rep(prior$mu.theta, H) 
    nu1[1] <- 1
    tau1[1] <- rgamma(1, prior$a, prior$b)
    th1[1] <- 0
    ga[1,]<- prior$mu.gamma 
  }
  else
  {
    w[1,] <- state$w
    nu0[1,] <- state$nu0
    nu1[1] <- 1
    th0[1,] <- state$th0
    tau0[1,] <- state$tau0
    th1[1] <- state$th1
    tau1[1] <- state$tau1
    ga[1,] <- state$ga
  }
  
  ## quantity of interest
  if(is.null(grid$grids))
  {
    x.grid <- seq(0, max.x, length=100)
    y.grid <- seq(min(y)-sqrt(var(y)), max(y)+sqrt(var(y)), length = 100)	
    beta_x <- matrix(NA, mcmc$nrep+mcmc$nb, length(x.grid))	
    f0 = matrix(NA, mcmc$nrep+mcmc$nb, length(y.grid))
    f1 = matrix(NA, mcmc$nrep+mcmc$nb, length(y.grid))
  }
  else
  {
    x.grid <- grid$xgrid
    y.grid <- grid$ygrid
    beta_x <- matrix(NA, mcmc$nrep+mcmc$nb, length(x.grid))	
    f0 = matrix(NA, mcmc$nrep+mcmc$nb, length(y.grid))
    f1 = matrix(NA, mcmc$nrep+mcmc$nb, length(y.grid))
  }
  
  z.val <- matrix(rep(z.val, length(y.grid)), ncol=length(z.val), byrow=T)
  
  ## basis expansion
  knots <- seq(min(x)+1, max.x, length=prior$J-3)
  basisX <- function(x) iSpline(x, df=3, knots = knots, Boundary.knots=c(0,max.x+1), intercept = FALSE)
  
  phiX <- basisX(x)
  phi.grid <- basisX(x.grid)
  
  ## beta_i is the interpolating function evaluated at x_i
  beta_i <- as.double(phiX %*% w[1, ])
  
  ## f0i and f1i
  f0i <- sapply(1:n, mixdensity_multi, y=y, z=z, nu=nu0[1,], theta=th0[1,], tau=tau0[1,], ga=ga[1,])
  f1i <- dnorm(y, rep(th1[1],n)+ crossprod(t(z),ga[1,]) , sqrt(1/tau1[1]))
  
  # start the MCMC simulation 
  set.seed(seed)
  for(ite in 2:(mcmc$nrep+mcmc$nb))
  {
    # 0. Print the iteration
    if(verbose)
    {
      if(ite==mcmc$nb) cat("Burn in done\n")
      if(ite %in% print_now) cat(ite, "iterations over",
                                 mcmc$nrep+mcmc$nb, "\n")
    }
    
    # 1. Update d_i marginalising out b_i from
    d <- rbinom(n, 1, prob=(beta_i*f1i)/((1-beta_i)*f0i + beta_i*f1i) )
    
    # 2. Update b_i from the multinomial 
    b <- sapply(1:n, labelling_b_multi, w[ite-1,], phi=phiX, f0i=f0i, f1i=f1i)
    
    # 3. Update c_i, marginalizing over b_i and d_i, from the multinomial 
    ind0 <- c(1:n)[d==0] #contiene l'indice delle oss nella componente a dose 0
    ind1 <- c(1:n)[d==1] #contiene l'indice delle oss nella componente a dose inf
    c <- sapply(ind0, labelling_c_multi, y=y, z=z, nu=nu0[ite-1,], 
                theta=th0[ite-1,], tau=tau0[ite-1,], ga=ga[ite-1,])
    
    # 4. Update the mixture weights sampling from dirichlet \per le oss t.c. d=0
    n_h <- table(factor(c, levels=1:H))
    nu0[ite,] <- as.double(rdirichlet(1, as.double(prior$alpha + n_h)))
    
    # 5. Update w from the Dirichlet and obtain an update function beta_i
    n_j <- table(factor(b, levels=1:J))  
    w[ite, ] <- as.double(rdirichlet(1, as.double(prior$eta + n_j)))
    
    beta_i <- as.numeric(phiX %*% w[ite, ])
    beta_i[beta_i>1] <- 1
    beta_i[beta_i<0] <- 0
    
    # 6. Update gamma
    
    ## in cluster 0: m=0
    indh <- sapply(c(1:H), function(x) c(1:length(ind0))[c==x]) # indh: lista di indici degli elementi di ind0 del gruppo h
    ZTZ <- lapply(indh, function(x) crossprod(z[ind0[x],,drop=F]) )
    # quantity 1: tau_0h*ZTZ for each h in group 0
    qt1 <- lapply(c(1:H), function(x) ZTZ[[x]]*tau0[ite-1,x])
    YTZ <- lapply(indh, function(x)  crossprod(y[ind0[x]], z[ind0[x],,drop=F]) )
    ITZ <- lapply(indh, function(x)  crossprod(rep(1,length(ind0[x])), z[ind0[x],,drop=F]) )
    # quantity 2: tau_0h*(YTZ-th_0h*Z) for each h in group 0
    qt2 <- lapply(c(1:H), function(x) tau0[ite-1,x]*(YTZ[[x]]- th0[ite-1,x]*ITZ[[x]]) )
    
    ## in cluster 1: m=1
    qt1_1 <- tau1[ite-1]*crossprod((z[ind1,,drop=F]))
    qt2_1 <- tau1[ite-1]*(crossprod(y[ind1],z[ind1,,drop=F]) - crossprod(rep(th1[ite-1],length(ind1)), z[ind1,,drop=F]))
    
    ## list of quantities in all possible clusters
    qt1[[H+1]] <- qt1_1
    qt2[[H+1]] <- qt2_1
    
    post.sigma <- solve( solve(prior$sigma.gamma)+Reduce('+', qt1) )
    
    post.mu <- tcrossprod( post.sigma, (crossprod(prior$mu.gamma, solve(prior$sigma.gamma)) + Reduce('+', qt2)) )
    
    ga[ite,] <- rmvnorm(1, mean=post.mu, sigma=post.sigma)

    # 7. Update theta and tau 
    n_0h <- n_h               
    n_1h <- length(d[d==1])
    
    ## in cluster 0: m=0
    hat_a <- prior$a + n_0h/2
    # quantity 3: (y-(th+Z ga))T(y-(th+Z ga))
    qt3 <- sapply(1:H, function(x) 
      crossprod( y[ind0[indh[[x]]]]-(rep(th0[ite-1,x],length(indh[[x]]))+crossprod(t(z[ind0[indh[[x]]],,drop=F]),ga[ite,])) ))
    hat_b <- prior$b + 0.5*(unlist(qt3))
    tau0[ite, ] <- rgamma(H, hat_a, hat_b)
    
    hat_kappa <- (1/prior$k.theta + n_0h*tau0[ite,])^-1
    qt4 <- lapply(indh, function(x) sum( t(y[ind0[x]]) - crossprod(ga[ite,], t(z[ind0[x],,drop=F])) ) )
    qt5 <- lapply(c(1:H), function(x) qt4[[x]]*tau0[ite,x])
    hat_mu <- hat_kappa * (prior$mu.theta/prior$k.theta + unlist(qt5) )
    th0[ite, ] <- rtruncnorm(H, a=th1[ite-1], b=Inf, hat_mu, sqrt(hat_kappa))
    
    ## in cluster 1: m=1
    hat_a <- prior$a + n_1h/2
    qt3 <- crossprod( ( y[ind1]- (rep(th1[ite-1],length(ind1)) + crossprod(t(z[ind1,,drop=F]),ga[ite,])) ) )
    hat_b <- prior$b + 0.5*(qt3)
    tau1[ite] <- rgamma(1, hat_a, hat_b)
    
    hat_kappa <- (1/prior$k.theta + n_1h*tau1[ite])^-1
    qt4 <- tau1[ite]*sum( t(y[ind1])-crossprod( ga[ite,],t(z[ind1,,drop=F]) ) )
    hat_mu <- hat_kappa * (prior$mu.theta/prior$k.theta + qt4)
    th1[ite] <- rtruncnorm(1, a=-Inf, b=min(th0[ite, ]), hat_mu, sqrt(hat_kappa))
    
    # update the values of the densities in the observed points
    f0i <- sapply(1:n, mixdensity_multi, y=y, z=z, nu=nu0[ite,], theta=th0[ite,], 
                  tau=tau0[ite,], ga=ga[ite,])
    f1i <- dnorm(y, (rep(th1[ite],n)+crossprod(t(z),ga[ite,])), sqrt(1/tau1[ite]))
    
    # 7. compute some posterior quantities of interest
    beta_x[ite, ] <- phi.grid %*% w[ite, ]
    f0[ite,] <- sapply(1:length(y.grid), mixdensity_multi, y=y.grid, z=z.val, 
                       nu=nu0[ite,], theta=th0[ite,], 
                       tau=tau0[ite,], ga=ga[ite,])
    f1[ite,] <- dnorm(y.grid, (th1[ite]+z.val%*%ga[ite,]) , sqrt(1/tau1[ite]))
    
  }
  
  post.mean.beta <- colMeans(beta_x[-c(1:mcmc$nb),])
  post.mean.w <- colMeans(w[-c(1:mcmc$nb),])
  post.mean.th0 <- colMeans(th0[-c(1:mcmc$nb),])
  post.mean.tau0 <- colMeans(tau0[-c(1:mcmc$nb),])
  post.mean.th1 <- mean(th1[-c(1:mcmc$nb)])
  post.mean.tau1 <- mean(tau1[-c(1:mcmc$nb)])
  post.mean.f0 <- colMeans(f0[-c(1:mcmc$nb),])
  post.mean.f1 <- colMeans(f1[-c(1:mcmc$nb),])
  post.mean.nu0 <- colMeans(nu0[-c(1:mcmc$nb),])
  post.mean.nu1 <- mean(nu1[-c(1:mcmc$nb)])
  post.mean.ga <- colMeans(ga[-c(1:mcmc$nb),,drop=F])

  ci.beta <- apply(beta_x[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
  ci.w<- apply(w[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
  ci.th0 <- apply(th0[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
  ci.tau0 <- apply(tau0[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
  ci.th1 <- quantile(th1[-c(1:mcmc$nb)], probs=c(0.025, 0.975))
  ci.tau1 <- quantile(tau1[-c(1:mcmc$nb)], probs=c(0.025, 0.975)) 
  ci.f0 <- apply(f0[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
  ci.f1 <- apply(f1[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
  ci.nu0 <- apply(nu0[-c(1:mcmc$nb),], 2, quantile, probs=c(0.025, 0.975))
  ci.ga <- apply(ga[-c(1:mcmc$nb),,drop=F], 2, quantile, probs=c(0.025, 0.975))

  # output
  output <- list(
    post.means = list(beta=post.mean.beta, w=post.mean.w, 
                     th0=post.mean.th0, tau0=post.mean.tau0, 
                     th1=post.mean.th1, tau1=post.mean.tau1, 
                     f0=post.mean.f0, f1=post.mean.f1,
                     nu0=post.mean.nu0, nu1=post.mean.nu1,
                     ga=post.mean.ga),
    ci = list( beta=ci.beta, w=ci.w,  
               th0=ci.th0, tau0=ci.tau0, 
               th1=ci.th1, tau1=ci.tau1, 
               f0=ci.f0, f1=ci.f1,
               nu0=ci.nu0, #nu1=ci.nu1,
               ga= ci.ga),
    mcmc = list(beta=beta_x, w=w, th0=th0, tau0=tau0, th1=th1, tau1=tau1,
                f0=f0, f1=f1, nu0=nu0, nu1=nu1, ga=ga)
    )
  list(out = output, z.val = z.val)
}
