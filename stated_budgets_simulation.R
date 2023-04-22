##################################################################################################
#####Simulation ##################################################################################
#####Including stated budgets (with measurement error) in the likelihood #########################
##################################################################################################

rm(list=ls())

# load packages
library(xtable)
library(devtools)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
library(bayesm)
library(ggplot2)
library(tikzDevice)
library(plyr)
library(latex2exp)

set.seed(5689484)

# individual information
cmax = 5 

# sample size
nunits = 1000
# variance-covariance matrix of betastars
sigma = diag(c(0.3,0.3,1,1),nrow=4,ncol=4)   
tsigma = t(sigma)
sigma[lower.tri(sigma)] = tsigma[lower.tri(tsigma)]

# mean values
avgbeta = c(1,0.5,3,1) 
# draw
betastar = mvrnorm(n=nunits, avgbeta, sigma) #draw true betas from assumed distribution of heterogeneity
# budget
beta = exp(betastar[,1])
# price as in BLP
beta = cbind(beta, exp(betastar[,2]))
# brands
beta = cbind(beta,betastar[,3:4])

# draw stated budgets (centered at true budgets with some error)
true_sd_error = 3
b_stated = array(0,dim=c(nunits,1))
for(n in 1:nunits){
  b_stated[n] = beta[n,1] + rnorm(1,mean=0,sd=true_sd_error)
}


### fraction of budget constrained individuals (max. p=4.5)
length(beta[beta[,1]<4.5,1])/nrow(beta)

###
# generate MNL data with BC
###
simmnl_BC_BLP = function(n,beta,budget,grid_1,grid_2){
  
  k = length(beta)
  # full matrix of product alts
  X_full_original = rbind(c(1,0),c(0,1),c(0,0))
  # deduce p
  p = dim(X_full_original)[1]
  # replicate original X-matrix
  for(i in 1:n){
    # sample prices from grid
    ind_1 = sample(seq(1,length(grid_1)),size=1)
    ind_2 = sample(seq(1,length(grid_2)),size=1)
    if(i==1){
      X_full = cbind(c(grid_1[ind_1],grid_2[ind_2],0),X_full_original)
    }else{
      X_full = rbind(X_full,cbind(c(grid_1[ind_1],grid_2[ind_2],0),X_full_original))
    }
  }
  # draw choices
  y=vector("double",n)
  ind=1:p
  # loop over tasks
  for(i in 1:n){
    constraint_ind = rep(1,p)
    X_BLP = X_full[((i-1)*p+1):(i*p),]
    for(np in 1:p){
      if(X_BLP[np,1] >= budget){
        constraint_ind[np] = 0 #set indicator to zero if price is larger than budget
      }else{
        X_BLP[np,1] = log(budget-X_BLP[np,1])
      }
    }
    Xbeta=X_BLP%*%beta
    Prob=exp(Xbeta)
    numerator = Prob*constraint_ind
    denominator = sum(numerator)
    Prob = numerator/denominator
    yvec=rmultinom(1,1,Prob)
    y[i]=ind%*%yvec
  }
  return(list(y=y,X=X_full,beta=beta,prob=Prob))
}

grid_1 = c(0.5,1,1.5,2,2.5,3,3.5,4,4.5) 
grid_2 = grid_1
lgtdata_sim=NULL
hdata=NULL
for (i in 1:nunits) {
  hdata[[i]] = simmnl_BC_BLP(cmax,beta[i,-1],beta[i,1],grid_1,grid_2)
}

# structure in a list
for(i in 1:nunits){
  lgtdata_sim[[i]]=list(X=hdata[[i]]$X,y=hdata[[i]]$y) 
}

E_Data = list(p=dim(lgtdata_sim[[1]]$X)[1]/cmax,lgtdata = lgtdata_sim)

#####################################
### model-free aggregate demand curve
#####################################

# price chosen per taks and exposure to max price in the experiment
p = nrow(E_Data$lgtdata[[1]]$X)/cmax
pr = 1
N = nunits

price_chosen_perTask <- array(0,dim=c(N,length(E_Data$lgtdata[[1]]$y)))
Exposed_max_price_overTasks <- array(0,dim=c(N,length(E_Data$lgtdata[[1]]$y)))
for(indi in 1:N){
  for(task in 1:length(E_Data$lgtdata[[1]]$y)){
    price_chosen_perTask[indi,task] = E_Data$lgtdata[[indi]]$X[((task-1)*p+E_Data$lgtdata[[indi]]$y[task]),pr]
    Exposed_max_price_overTasks[indi,task] = if(max(E_Data$lgtdata[[indi]]$X[((task-1)*p+1):(task*p),pr])==max(grid_1)){1}else{0}
  }
}
max_price_chosen = apply(price_chosen_perTask,1,max)
Exposed_max_price_overTasks = apply(Exposed_max_price_overTasks,1,sum)

summary(max_price_chosen)
summary(Exposed_max_price_overTasks)

##############################################################
#######Run BLP Budget model using stated budget###############
##############################################################
Rcpp::sourceCpp("functions/main_mcmc_loop_stated_budgets.cpp",showOutput = FALSE)
source('functions/main_function_stated_budgets.R')

# number of constrained coefficients (budget & price)
nvar_c = 2
# position of price coefficient in design matrix
pr=1

# prior setting (Pachali, Kurz and Otter, QME 2020)
Amu = diag(1/10, nrow = nvar_c, ncol = nvar_c)
mustarbarc = matrix(rep(0, nvar_c), nrow = nvar_c)
nu = 15 + nvar_c
V = nu * diag(nvar_c)*0.5

Prior = list(ncomp=1, Amu = Amu, mustarbarc = mustarbarc, nu = nu, V = V)
Mcmc = list(R=10000, keep=1)
max_p = max(grid_1)
true_betastar = betastar
out_BC = rhierMnlRwMixture_SR(Data=E_Data,stated_budget=b_stated,Prior=Prior,Mcmc=Mcmc,nvar_c=nvar_c,pr=pr,starting_budget = log(max_p+0.1))

betastar_HB_BC = out_BC$betadraw
compdraw_HB = out_BC$nmix$compdraw
probdraw_HB = out_BC$nmix$probdraw
rejection = out_BC$rejection
loglike_BC = out_BC$loglike
error_sd = out_BC$error_sd

# compute rejection rate of sampler 
rej_rate_indi = apply(rejection,2,mean)
summary(rej_rate_indi)
rej_rate_agg = mean(rej_rate_indi)

R=dim(betastar_HB_BC)[3]

burnin = 8000
R = dim(betastar_HB_BC)[3]

betastar_HB_BC = betastar_HB_BC[,,(burnin+1):R]
compdraw_HB = compdraw_HB[(burnin+1):R]
probdraw_HB = probdraw_HB[(burnin+1):R]
rejection = rejection[(burnin+1):R,]
loglike_BC = loglike_BC[(burnin+1):R]
error_sd = error_sd[(burnin+1):R]

R = dim(betastar_HB_BC)[3]
N = dim(betastar_HB_BC)[1]

### recovery of dgp
mean(error_sd)
quantile(error_sd,probs = c(0.025,0.975))
true_sd_error

par(mfrow=c(1,1))  # multiple plots are filled by rows!!!
plot(error_sd,xlab="R",ylab="sd error", ylim=c(0.5,5) ,type="l",cex.main=2,cex.lab=1.5 ,font.lab=2 ,font.axis=1 ,cex.axis=1.5);grid()
abline(h=true_sd_error,col="red",lwd=5)

### recovery of beta
# estimated
betastar_est <- array(aperm(betastar_HB_BC, perm=c( 1 , 3 , 2 )),
                      dim=c(dim(betastar_HB_BC)[1] * dim(betastar_HB_BC)[3],
                            dim(betastar_HB_BC)[2]))
beta_est = betastar_est
beta_est[,1] = exp(betastar_est[,1]) # budget
beta_est[,2] = exp(betastar_est[,2]) # BLP-price
summary(beta_est)
# true
summary(beta)

### check error (RMSE) b/w posterior mean and true budgets at the individual level
beta_budgets_draws = exp(betastar_HB_BC[,1,])
beta_budgets_PM = apply(beta_budgets_draws,1,mean)
beta_budgets_SD = apply(beta_budgets_draws,1,sd)

RMSE_stated_budgets = sqrt(mean((beta_budgets_PM-beta[,1])^2))

RMSE_stated_budgets
mean(beta_budgets_SD)

#################################################
### Infer budgets w/o using stated budgets ######
#################################################

# load budget sampler now...
Rcpp::sourceCpp("functions/w_o_stated_budgets/rhierMnlRwMixture_rcpp_loop_Illus_BLP_type.cpp",showOutput = FALSE)
source('functions/w_o_stated_budgets/rhierMnlRwMixture_main_BC.R')

# number of constrained coefficients (budget & price)
nvar_c = 2
# position of price coefficient in design matrix
pr=1

# prior setting 
Amu = diag(1/10, nrow = nvar_c, ncol = nvar_c)
mustarbarc = matrix(rep(0, nvar_c), nrow = nvar_c)
nu = 15 + nvar_c
V = nu * diag(nvar_c)*0.5

Prior = list(ncomp=1, Amu = Amu, mustarbarc = mustarbarc, nu = nu, V = V)
out_BC = rhierMnlRwMixture_SR(Data=E_Data,Prior=Prior,Mcmc=Mcmc,nvar_c=nvar_c,pr=pr,starting_budget = log(max_p+0.1))

betastar_HB_BC = out_BC$betadraw
compdraw_HB = out_BC$nmix$compdraw
probdraw_HB = out_BC$nmix$probdraw
rejection = out_BC$rejection
loglike_BC = out_BC$loglike

# compute rejection rate of sampler 
rej_rate_indi = apply(rejection,2,mean)
summary(rej_rate_indi)
rej_rate_agg = mean(rej_rate_indi)

R=dim(betastar_HB_BC)[3]

burnin = 8000
R = dim(betastar_HB_BC)[3]

betastar_HB_BC = betastar_HB_BC[,,(burnin+1):R]

### check error (RMSE) b/w posterior mean and true budgets at the individual level
beta_budgets_draws_only = exp(betastar_HB_BC[,1,])
beta_budgets_PM_only = apply(beta_budgets_draws_only,1,mean)
beta_budgets_SD_only = apply(beta_budgets_draws_only,1,sd)

RMSE_only_budgets = sqrt(mean((beta_budgets_PM_only-beta[,1])^2))

# RMSE only estimating budgets
RMSE_only_budgets
# RMSE including stated budgets
RMSE_stated_budgets

# mean posterior standard deviations (over individuals) only estimating budgets
mean(beta_budgets_SD_only)
# mean posterior standard deviations (over individuals) including stated budgets
mean(beta_budgets_SD)














