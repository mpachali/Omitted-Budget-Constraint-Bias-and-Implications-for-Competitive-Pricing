rhierMnlRwMixture_SR <- function (Data, Prior, Mcmc, nvar_c, pr, starting_budget) 
{
  fsh=function() 
  {
    # 
    # P. Rossi
    # revision history: 3/27/05
    #
    # Purpose:
    #  function to flush console (needed only under windows)
    #
    if (Sys.info()[1] == "Windows") flush.console()
    return()
  }
  
  if (missing(Data)) {
    pandterm("Requires Data argument -- list of p,lgtdata, and (possibly) Z")
  }
  if (is.null(Data$p)) {
    pandterm("Requires Data element p (# chce alternatives)")
  }
  p = Data$p
  if (is.null(Data$lgtdata)) {
    pandterm("Requires Data element lgtdata (list of data for each unit)")
  }
  lgtdata = Data$lgtdata
  nlgt = length(lgtdata)
  drawdelta = TRUE
  if (is.null(Data$Z)) {
    cat("Z not specified", fill = TRUE)
    fsh()
    drawdelta = FALSE
  }
  else {
    if (nrow(Data$Z) != nlgt) {
      pandterm(paste("Nrow(Z) ", nrow(Z), "ne number logits ", 
                     nlgt))
    }
    else {
      Z = Data$Z
    }
  }
  if (drawdelta) {
    nz = ncol(Z)
    colmeans = apply(Z, 2, mean)
    if (sum(colmeans) > 1e-05) {
      pandterm(paste("Z does not appear to be de-meaned: colmeans= ", 
                     colmeans))
    }
  }
  ypooled = NULL
  Xpooled = NULL
  if (!is.null(lgtdata[[1]]$X)) {
    oldncol = ncol(lgtdata[[1]]$X) + 1
  }
  for (i in 1:nlgt) {
    if (is.null(lgtdata[[i]]$y)) {
      pandterm(paste("Requires element y of lgtdata[[", 
                     i, "]]"))
    }
    if (is.null(lgtdata[[i]]$X)) {
      pandterm(paste("Requires element X of lgtdata[[", 
                     i, "]]"))
    }
    ypooled = c(ypooled, lgtdata[[i]]$y)
    nrowX = nrow(lgtdata[[i]]$X)
    if ((nrowX/p) != length(lgtdata[[i]]$y)) {
      pandterm(paste("nrow(X) ne p*length(yi); exception at unit", 
                     i))
    }
    newncol = ncol(lgtdata[[i]]$X) + 1 
    if (newncol != oldncol) {
      pandterm(paste("All X elements must have same # of cols; exception at unit", 
                     i))
    }
    Xpooled = rbind(Xpooled, lgtdata[[i]]$X)
    oldncol = newncol
  }
  nvar = ncol(Xpooled) + 1 ###Budget constraint parameter not in design
  nvar_uc = nvar-nvar_c
  levely = as.numeric(levels(as.factor(ypooled)))
  if (length(levely) != p) {
    pandterm(paste("y takes on ", length(levely), " values -- must be = p"))
  }
  bady = FALSE
  for (i in 1:p) {
    if (levely[i] != i) 
      bady = TRUE
  }
  cat("Table of Y values pooled over all units", fill = TRUE)
  print(table(ypooled))
  if (bady) {
    pandterm("Invalid Y")
  }
  if (missing(Prior)) {
    pandterm("Requires Prior list argument (at least ncomp)")
  }
  if (is.null(Prior$ncomp)) {
    pandterm("Requires Prior element ncomp (num of mixture components)")
  }
  else {
    ncomp = Prior$ncomp
  }
  if (is.null(Prior$mustarbarc)) {
    mustarbarc = matrix(rep(0, nvar_c), nrow = nvar_c)
  }
  else {
    mustarbarc = matrix(Prior$mustarbarc, nrow = nvar_c)
  }
  if (is.null(Prior$Gammabar)) {
    Gammabar = matrix(0, nrow = (nvar_c+1), ncol = nvar_uc)
  }
  else {
    Gammabar = matrix(Prior$Gammabar, nrow = (nvar_c+1), ncol = nvar_uc)
  }
  if (is.null(Prior$Amu)) {
    Amu = diag(0.01, nrow = nvar_c, ncol = nvar_c)
  }
  else {
    Amu = matrix(Prior$Amu, nrow = nvar_c, ncol = nvar_c)
  }
  if (is.null(Prior$A_Gamma)) {
    A_Gamma = diag(0.01, nrow = (nvar_c+1), ncol = (nvar_c+1))
  }
  else {
    A_Gamma = matrix(Prior$A_Gamma, nrow = (nvar_c+1), ncol = (nvar_c+1))
  }
  if (is.null(Prior$nu)) {
    nu = nvar_c + 15
    ###Informative prior on restricted coefficients
  }
  else {
    nu = Prior$nu
  }
  if (is.null(Prior$nu_Sigma)) {
    nu_Sigma = nvar_uc + 5
  }
  else {
    nu_Sigma = Prior$nu_Sigma
  }
  if (is.null(Prior$V)) {
    V = nu * diag(nvar_c)*0.5
    ###Allenby et al. (2014) 
  }
  else {
    V = Prior$V
  }
  if (is.null(Prior$V_Sigma)) {
    V_Sigma = nu_Sigma * diag(nvar_uc)
  }
  else {
    V_Sigma = Prior$V_Sigma
  }
  if (is.null(Prior$Ad) & drawdelta) {
    Ad = 0.01 * diag(nvar * nz)
  }
  else {
    Ad = Prior$Ad
  }
  if (drawdelta) {
    if (ncol(Ad) != nvar * nz | nrow(Ad) != nvar * nz) {
      pandterm("Ad must be nvar*nz x nvar*nz")
    }
  }
  if (is.null(Prior$deltabar) & drawdelta) {
    deltabar = rep(0, nz * nvar)
  }
  else {
    deltabar = Prior$deltabar
  }
  if (drawdelta) {
    if (length(deltabar) != nz * nvar) {
      pandterm("deltabar must be of length nvar*nz")
    }
  }
  if (is.null(Prior$a)) {
    a = rep(5, ncomp)
  }
  else {
    a = Prior$a
  }
  if (length(a) != ncomp) {
    pandterm("Requires dim(a)= ncomp (no of components)")
  }
  bada = FALSE
  for (i in 1:ncomp) {
    if (a[i] < 0) 
      bada = TRUE
  }
  if (bada) 
    pandterm("invalid values in a vector")
  if (missing(Mcmc)) {
    pandterm("Requires Mcmc list argument")
  }
  else {
    if (is.null(Mcmc$s)) {
      s = 2.38/sqrt(nvar)
    }
    else {
      s = Mcmc$s
    }
    if (is.null(Mcmc$w)) {
      w = 0.1
    }
    else {
      w = Mcmc$w
    }
    if (is.null(Mcmc$keep)) {
      keep = 1
    }
    else {
      keep = Mcmc$keep
    }
    if (is.null(Mcmc$R)) {
      pandterm("Requires R argument in Mcmc list")
    }
    else {
      R = Mcmc$R
    }
    if (is.null(Mcmc$nprint)) {
      nprint = 100
    }
    else {
      nprint = Mcmc$nprint
    }
    if (nprint < 0) {
      pandterm("nprint must be an integer greater than or equal to 0")
    }
  }
  cat(" ", fill = TRUE)
  cat("Starting MCMC Inference for Hierarchical Logit:", fill = TRUE)
  cat("   Normal Mixture with", ncomp, "components for first stage prior", 
      fill = TRUE)
  cat(paste("  ", p, " alternatives; ", nvar, " variables in X"), 
      fill = TRUE)
  cat(paste("   for ", nlgt, " cross-sectional units"), fill = TRUE)
  cat(" ", fill = TRUE)
  cat("Prior Parms: ", fill = TRUE)
  cat("nu =", nu, "nu_Sigma =", nu_Sigma,  fill = TRUE)
  cat("V ", fill = TRUE)
  print(V)
  cat("V_Sigma ", fill = TRUE)
  print(V_Sigma)
  cat("mustarbarc ", fill = TRUE)
  print(mustarbarc)
  cat("Amu ", fill = TRUE)
  print(Amu)
  cat("Gammabar ", fill = TRUE)
  print(Gammabar)
  cat("A_Gamma ", fill = TRUE)
  print(A_Gamma)
  cat("a ", fill = TRUE)
  print(a)
  if (drawdelta) {
    cat("deltabar", fill = TRUE)
    print(deltabar)
    cat("Ad", fill = TRUE)
    print(Ad)
  }
  cat(" ", fill = TRUE)
  cat("MCMC Parms: ", fill = TRUE)
  cat(paste("s=", round(s, 3), " w= ", w, " R= ", R, " keep= ", 
            keep, " nprint= ", nprint), fill = TRUE)
  cat("", fill = TRUE)
  oldbetas = matrix(double(nlgt * nvar), ncol = nvar)
  ind_stay = array(0,dim=c((R/keep),nlgt))
  ind_stay_temp = array(0,dim=c(nlgt,1))
  
  llmnlFract = function(beta, y, X, betapooled, rootH, w, wgt, pr) {
    z = as.vector(rootH %*% (beta - betapooled))
    return((1 - w) * llmnl_BC(beta, y, X, pr) + w * wgt * (-0.5 * 
                                                             (z %*% z)))
  }
  
  ###################################################
  ### Needs to be adjusted for BLP-type utility model
  ###################################################
  
  mnlHess_BC = function (betastar, y, X, pr) {
    beta = startobeta(betastar)  # call to compiled transformation routine
    n = length(y)
    j = nrow(X)/n
    k = ncol(X)
    beta_part = beta[-1]
    budget = beta[1]
    Hess = matrix(double(k * k), ncol = k)
    for (i in 1:n) {
      constraint_ind = rep(1,j)
      for(np in 1:j){
        if(X[((i-1)*j+np),pr] >= budget){
          constraint_ind[np] = 0 #set indicator to zero if price is larger than budget
        }
      }
      Xbeta=X[((i-1)*j+1):(i*j),]%*%beta_part
      ###Stabilization: out-of-budget & Xbeta --> infti at same time
      Xbeta=Xbeta*constraint_ind
      Prob=exp(Xbeta)
      numerator = Prob*constraint_ind
      denominator = sum(numerator)
      Prob = numerator/denominator
      Prob = as.vector(Prob)
      #p = as.vector(Prob[i, ])
      #A = diag(p) - outer(p, p)
      A = diag(Prob) - outer(Prob, Prob)
      Xt = X[(j * (i - 1) + 1):(j * i), ]
      Hess = Hess + crossprod(Xt, A) %*% Xt
    }
    Hess_full = matrix(0,k+1,ncol=k+1)
    Hess_full[1,1] = 1
    Hess_full[-1,-1] = Hess
    return(Hess_full)
  }

  ##################################################################################################
  ###TUNE CANDIDATES USING THE APPROXIMATIVE HESSIAN FOR BETASTARS##################################
  ##################################################################################################
  cat("initializing Metropolis candidate densities for ", nlgt, 
      " units ...", fill = TRUE)
  fsh()
  ###Initialize betastar 
  betastarinit = c(starting_budget,rep(0, nvar-1))
  ###Obtain Mle of betastar from pooled data
  # out = optim(betastarinit, llmnl_BC, method = "Nelder-Mead", control = list(fnscale = -1, 
  #                                                                            trace = 0, reltol = 1e-6, maxit = 200000), X = Xpooled, y = ypooled, pr = pr)
  # if(out$convergence != 0){
  #   pandterm("Optim has not converged!")
  # }else{
  #   betastarpooled = out$par
  # }
  
  betastarpooled = betastarinit
  
  ############################################
  ###Hessian of betastarpooled################
  ############################################
  # HESS_betapooled = mnlHess_BC(betastarpooled,ypooled,Xpooled,pr)
  # require(pracma)
  # HESS_betastarpooled = t(jacobian(startobeta, betastarpooled))%*%HESS_betapooled%*%jacobian(startobeta, betastarpooled)
  # rootHESS_betastarpooled = chol(HESS_betapooled)%*%jacobian(startobeta, betastarpooled) ###allows for more stable computation of root
  
  ############################################
  ##Individual hessians#######################
  ############################################
  
  for (i in 1:nlgt) {
    # wgt = length(lgtdata[[i]]$y)/length(ypooled)
    # out = optim(betastarpooled, llmnlFract, method = "Nelder-Mead", 
    #             control = list(fnscale = -1, trace = 0, reltol = 1e-5, maxit = 200000), 
    #             X = lgtdata[[i]]$X, y = lgtdata[[i]]$y, pr = pr, betapooled = betastarpooled, 
    #             rootH = rootHESS_betastarpooled, w = w, wgt = wgt, hessian = FALSE)
    # ###if convergence is achieved...
    # if (out$convergence == 0) {
    #   betafmle_temp = out$par
    #   if(betafmle_temp[1]>starting_budget){
    #     betafmle_temp[1] = starting_budget + 0.1*starting_budget
    #   }
    #   hess = mnlHess_BC(betafmle_temp, lgtdata[[i]]$y, lgtdata[[i]]$X, pr)
    #   hess_star = t(jacobian(startobeta, betafmle_temp))%*%hess%*%jacobian(startobeta, betafmle_temp)
    #   lgtdata[[i]] = c(lgtdata[[i]], list(converge = 1, 
    #                                       betafmle = betafmle_temp, hess = hess_star))
    # }
    # else {
      lgtdata[[i]] = c(lgtdata[[i]], list(converge = 0, 
                                          betafmle = c(starting_budget,rep(0, nvar-1)), hess = 1*diag(nvar)))
    # }
    oldbetas[i, ] = lgtdata[[i]]$betafmle
    if (i%%50 == 0) 
      cat("  completed unit #", i, fill = TRUE)
    fsh()
  }
  
  ###################
  ###END OF TUNING###
  ###################
  
  ind = NULL
  ninc = floor(nlgt/ncomp)
  for (i in 1:(ncomp - 1)) {
    ind = c(ind, rep(i, ninc))
  }
  if (ncomp != 1) {
    ###if there is more than a single component, equally divide individuals across components
    ind = c(ind, rep(ncomp, nlgt - length(ind)))
  }
  else {
    ind = rep(1, nlgt)
  }
  oldprob = rep(1/ncomp, ncomp)
  if (drawdelta) {
    olddelta = rep(0, nz * nvar)
  }
  else {
    olddelta = 0
    Z = matrix(0)
    deltabar = 0
    Ad = matrix(0)
  }
  draws = rhierMnlRwMixture_rcpp_loop_mod(lgtdata, nvar_c, Z, deltabar, Ad, mustarbarc, Amu, nu, V, Gammabar, A_Gamma, nu_Sigma, V_Sigma,
                                          s, R, keep, nprint, drawdelta, as.matrix(olddelta), a, oldprob, oldbetas, ind, ind_stay_temp, ind_stay,ncomp, pr)
  # draws = rhierMnlRwMixture_rcpp_loop_mod(lgtdata, nvar_c, Z, deltabar, Ad, mustarbarc, Amu, nu, V, Gammabar, A_Gamma, nu_Sigma, V_Sigma,
  #                                         s, R, keep, nprint, drawdelta, as.matrix(olddelta), a, oldprob, oldbetas, ind, betas_tuning, 
  #                                         betastarpooled, ind_stay_temp, ind_stay,ncomp,pr)
  
  if (drawdelta) {
    attributes(draws$Deltadraw)$class = c("bayesm.mat", "mcmc")
    attributes(draws$Deltadraw)$mcpar = c(1, R, keep)
  }
  attributes(draws$betadraw)$class = c("bayesm.hcoef")
  attributes(draws$nmix)$class = "bayesm.nmix"
  return(draws)
}
