rm(list=ls()) #Clean-up working environment

##################################################
# Install (if needed) and load required packages #
##################################################

if(!requireNamespace('metafor', quietly=TRUE)) { install.packages('metafor') }; library(lmtest)
if(!requireNamespace('lmtest', quietly=TRUE)) { install.packages('lmtest') }; library(lmtest)
if(!requireNamespace('doParallel', quietly=TRUE)) { install.packages('doParallel') }; library(doParallel)

##########################
# User-defined functions #
##########################

#MR-Egger regression
MREgger <- function(BetaXG,        #BetaXG:      SNP-exposure coefficient, typically a linear (if the exposure is continuous) or logistic (if the exposure is binary) regression coefficient.
                    BetaYG,        #BetaYG:      SNP-outcome coefficient, typically a linear (if the outcome is continuous) or logistic (if the outcome is binary) regression coefficient.
                    seBetaXG,      #seBetaXG:    standard error of BetaXG.
                    seBetaYG,      #seBetaYG:    standard error of BetaYG.
                    weighting,     #weighting:   unweighted (0), 1st order weights (1), modified weights (2).
                    alpha=0.05,    #alpha:       1-alpha is the confidence level of confidence intervals.
                    recode,        #recode:      if TRUE, ensure all estimated SNP-X coefficients are positive; otherwise (i.e., if FALSE), do nothing
                    remove_weak) { #remove_weak: if TRUE, remove the 10% of the variants with the smallest abs(BetaXG/seBetaXG) to mitigate incorrect positive coding; otherwise (i.e., if FALSE), do nothing.

  if(remove_weak) {
    ind_remove <- abs(BetaXG/seBetaXG)<sort(abs(BetaXG/seBetaXG))[round(length(BetaXG)*0.1)]
    
    BetaXG   <- BetaXG[!ind_remove]
    BetaYG   <- BetaYG[!ind_remove]
    seBetaXG <- seBetaXG[!ind_remove]
    seBetaYG <- seBetaYG[!ind_remove]
  }
  
  #Calculate weights
  if(weighting==0) {
    weights <- rep(1,length(seBetaYG)) #Unweighted
    
  } else if(weighting%in%c(1,2)) {
    weights <- 1/seBetaYG^2 #1st order weights
    
  } else {
    break('weighting must be 0, 1 or 2')
  }

  if(recode) {
    #Pre-processing steps to ensure all instrument-exposure estimates are positive
    BYG <- BetaYG*sign(BetaXG)
    BXG <- abs(BetaXG)  
  } else {
    BYG <- BetaYG
    BXG <- BetaXG
  }

  #MR-Egger model
  MREggerFit <- summary(lm(BYG ~ BXG, weights=weights))
  
  if(weighting==2) {
    weights    <- 1/(seBetaYG^2 + MREggerFit$coef[2,1]^2*seBetaXG^2)
    MREggerFit <- summary(lm(BYG ~ BXG, weights=weights))
  }
  
  #Inference with correct standard errors
  MREggerBeta0 <- MREggerFit$coef[1,1]
  MREggerBeta1 <- MREggerFit$coef[2,1]
  MREggerSE0   <- MREggerFit$coef[1,2]
  MREggerSE1   <- MREggerFit$coef[2,2]
  
  if(weighting>0) {
    MREggerSE0 <- MREggerSE0/min(MREggerFit$sigma, 1)
    MREggerSE1 <- MREggerSE1/min(MREggerFit$sigma, 1)
  }
  
  DF         <- length(BetaYG)-2
  MRBeta0_CI <- MREggerBeta0 + c(-1,1)*qt(df=DF, 1-alpha/2)*MREggerSE0
  MRBeta1_CI <- MREggerBeta1 + c(-1,1)*qt(df=DF, 1-alpha/2)*MREggerSE1
  MRBeta0_p  <- 2*(1-pt(abs(MREggerBeta0/MREggerSE0),DF))
  MRBeta1_p  <- 2*(1-pt(abs(MREggerBeta1/MREggerSE1),DF))

  MREggerResults <- matrix(nrow=2, ncol=5)
  MREggerResults[1,] <- c(MREggerBeta0, MREggerSE0, MRBeta0_CI ,MRBeta0_p)  
  MREggerResults[2,] <- c(MREggerBeta1, MREggerSE1, MRBeta1_CI, MRBeta1_p)
  colnames(MREggerResults) <- c('Estimate', 'SE', 'CI_low', 'CI_upp', 'P-value')
  rownames(MREggerResults) <- c('Intercept', 'Beta')  
  
  return(MREggerResults)
}

#Heterogeneity test
Het_test <- function(BetaXG, BetaYG, seBetaXG, seBetaYG, weighting, recode, remove_weak) {

  if(remove_weak) {
    ind_remove <- abs(BetaXG/seBetaXG)<sort(abs(BetaXG/seBetaXG))[round(length(BetaXG)*0.1)]
    
    BetaXG   <- BetaXG[!ind_remove]
    BetaYG   <- BetaYG[!ind_remove]
    seBetaXG <- seBetaXG[!ind_remove]
    seBetaYG <- seBetaYG[!ind_remove]
  }
  
  #Calculate weights
  if(weighting==0) {
    weights <- rep(1,length(seBetaYG)) #Unweighted
    
  } else if(weighting%in%c(1,2)) {
    weights <- 1/seBetaYG^2 #1st order weights
    
  } else {
    break('weighting must be 0, 1 or 2')
  }  
  
  #Cochran's Q
  IVWfit <- summary(lm(BetaYG ~ -1+BetaXG, weights=weights)) #IVW model
  Q_DF   <- length(BetaXG)-1                                 #Degrees of freedom
  
  if(weighting==2) {
    Q_stat <- sum((BetaYG - IVWfit$coef[1,1]*BetaXG)^2 / (seBetaYG^2 + IVWfit$coef[1,1]^2*seBetaXG^2)) #Q statistic
  } else {
    Q_stat <- Q_DF*IVWfit$sigma^2 #Q statistic
  }
  
  Q_I2 <- max(0, (Q_stat-Q_DF)/Q_stat)*100 #I2
  Q_p  <- 1-pchisq(Q_stat, Q_DF)           #P-value
  
  #Rucker's Q'
  
  if(recode) {
    #Pre-processing steps to ensure all instrument-exposure estimates are positive
    BYG <- BetaYG*sign(BetaXG)
    BXG <- abs(BetaXG)  
  } else {
    BYG <- BetaYG
    BXG <- BetaXG
  }
  
  MREggerFit <- summary(lm(BYG ~ BXG, weights=weights)) #MR-Egger model
  Qline_DF   <- length(BXG)-2                           #Degrees of freedom
  
  if(weighting==2) {
    Qline_stat <- sum((BYG - MREggerFit$coef[2,1]*BXG - MREggerFit$coef[1,1])^2 / (seBetaYG^2 + MREggerFit$coef[2,1]^2*seBetaXG^2)) #Q' statistic
  } else {
    Qline_stat <- Qline_DF*MREggerFit$sigma^2 #Q' statistic
  }
  
  Qline_I2 <- max(0, (Qline_stat-Qline_DF)/Qline_stat)*100 #I2
  Qline_p  <- 1-pchisq(Qline_stat, Qline_DF)               #P-value

  return(c(Q_stat,     Q_DF,     Q_I2,     Q_p,
           Qline_stat, Qline_DF, Qline_I2, Qline_p))
}

#Heteroscedasticity test
Var_test <- function(BetaXG, BetaYG, seBetaXG, seBetaYG, weighting, recode, remove_weak) {

  if(remove_weak) {
    ind_remove <- abs(BetaXG/seBetaXG)<sort(abs(BetaXG/seBetaXG))[round(length(BetaXG)*0.1)]
    
    BetaXG   <- BetaXG[!ind_remove]
    BetaYG   <- BetaYG[!ind_remove]
    seBetaXG <- seBetaXG[!ind_remove]
    seBetaYG <- seBetaYG[!ind_remove]
  }
  
  #Calculate weights
  if(weighting==0) {
    weights <- rep(1,length(seBetaYG)) #Unweighted
    
  } else if(weighting%in%c(1,2)) {
    weights <- 1/seBetaYG^2 #1st order weights
    
  } else {
    break('weighting must be either T or F')
  }
  
  if(recode) {
    #Pre-processing steps to ensure all instrument-exposure estimates are positive
    BYG <- BetaYG*sign(BetaXG)
    BXG <- abs(BetaXG)  
  } else {
    BYG <- BetaYG
    BXG <- BetaXG
  }
  
  weights <- length(weights)*weights/sum(weights) #Normalize weights so they sum to the number of genetic variants
  
  #Heteroscedasticity test
  fit <- lm(BYG ~ BXG, weights=weights) #MR-Egger model
  
  if(weighting==2) {
    weights <- 1/(seBetaYG^2 + fit$coef[2]^2*seBetaXG^2)
    weights <- length(weights)*weights/sum(weights)
    fit     <- summary(lm(BYG ~ BXG, weights=weights))
  }
  
  pval1 <- bptest(fit, varformula = ~ BXG)$p.value
  pval2 <- bptest(fit, varformula = ~ BXG + I(BXG^2))$p.value
  
  return(c(pval1,pval2))
}

#Function to extract desired quantities from a simulated dataset
get.res <- function(x, beta) {
  bias     <- x[1]-beta
  se       <- x[2]
  coverage <- x[3]<=beta & x[4]>=beta
  power    <- x[3]>0 | x[4]<0
  return(c(bias, se, coverage, power))
}

#########################
# Set up the simulation #
#########################

#Set up parallel computing
library(doParallel)
registerDoParallel(cores=n.cores) # Replace n.cores with the number of cores too be used in the simulations

n.sim <- 1e4 #Number of simulated datasets per scenario

#Parameters that will vary within scenarios
L   <- c(50,100,150,200)          #Number of candidate genetic instruments
n_X <- seq(2.5e4,2e5,2.5e4)       #Sample size where instrument-exposure associations were estimated
L   <- rep(L,each=length(n_X))    
n_X <- rep(n_X,length(unique(L)))
n_Y <- n_X                        #Sample size where instrument-outcome associations were estimated

#Define simulations scenarios
scenarios <- data.frame(scenario=1:6,                #Scenarios 1-2, 3 and 4-6 here respectively correspond to scenarios 1, 2 and 3 in the text
                        beta=c(0,0.5,0,0,0,0),       #Causal effect of the exposure on the outcome
                        rho_0=c(0,0,1,0,0,0),        #Proportion of variants with horizontal pleiotropic effects under the INSIDE assumption
                        rho_psi=c(0,0,0,0.5,0.75,1)) #Proportion of variants with horizontal pleiotropic effects violating the INSIDE assumption

scenarios <- scenarios[rep(1:nrow(scenarios), each=length(L)),]
scenarios$L <- rep(L,length(unique(scenarios$scenario)))
scenarios$n_X <- rep(n_X,length(unique(scenarios$scenario)))
scenarios$n_Y <- rep(n_Y,length(unique(scenarios$scenario)))

sim.res <- NULL

###################
# Run simulations #
###################

for(a in 1:nrow(scenarios)) {
  
  print(a)
  
  L        <- scenarios$L[a]
  n_X      <- scenarios$n_X[a]
  n_Y      <- scenarios$n_Y[a]
  scenario <- scenarios$scenario[a]
  beta     <- scenarios$beta[a]
  rho_0    <- scenarios$rho_0[a]
  rho_psi  <- scenarios$rho_psi[a]
  
  sim.res.a <- foreach(b=1:n.sim, .combine='rbind', .inorder=F, .packages=c('metafor','lmtest')) %dopar% {

    #Effect allele frequency (theta)
    theta  <- runif(L, 0.1, 0.9)
    maf    <- theta
    maf[theta>0.5] <- 1-theta[theta>0.5]
    
    min_delta <- 0.01
    max_delta <- 0.082
    delta     <- runif(L, min_delta, max_delta)                                    #Direct effect of each variant on X
    I_rho     <- sample(c(rep(0, round(L*(1-rho_psi))), round(rep(1, L*rho_psi))))
    phi       <- I_rho*runif(L, min_delta, max_delta)                              #Direct effect of each variant on U
    kappa_X   <- runif(L, 0.05, 0.5)
    kappa_Y   <- kappa_X
    
    I_0     <- sample(c(rep(0, round(L*(1-rho_0))), round(rep(1, L*rho_0))))
    alpha_0 <- I_0*runif(L, min_delta, max_delta) #Horizontal pleiotropy term (under INSIDE)
    
    #Variance (sigma2) of each genetic variant
    sigma2 <- 2*theta*(1-theta)

    gamma     <- delta + phi*kappa_X #Total effect of each variant on X

    lambda_X  <- gamma^2*sigma2                  #Amount of variance in X explained by each variant
    sigma_X   <- sqrt((1-lambda_X)/(sigma2*n_X)) #Standard error associated with gamma
    gamma_hat <- rnorm(L, gamma, sigma_X)        #Finite-sample estimate of gamma
    
    Gamma     <- alpha_0 + gamma*beta + phi*kappa_Y #Effect of each variant on Y
    lambda_Y  <- Gamma^2*sigma2                     #Amount of variance in Y explained by each genetic variant 
    sigma_Y   <- sqrt((1-lambda_Y)/(sigma2*n_Y))    #Standard error associated with Gamma
    Gamma_hat <- rnorm(L, Gamma, sigma_Y)           #Finite-sample estimate of Gamma

    #Run analyses and collect results

    #MR-Egger regression
    MREgger.res1 <- MREgger(gamma_hat, Gamma_hat, sigma_X, sigma_Y, weighting=1, recode=T, remove_weak=F)[2,] #First-order weights, estimated SNP-X coefficients, recode to all-positive
    MREgger.res2 <- MREgger(gamma_hat, Gamma_hat, sigma_X, sigma_Y, weighting=2, recode=T, remove_weak=F)[2,] #Modified weights, estimated SNP-X coefficients, recode to all-positive
    MREgger.res3 <- MREgger(gamma,     Gamma_hat, sigma_X, sigma_Y, weighting=1, recode=T, remove_weak=F)[2,] #First-order weights, true SNP-X coefficients, recode to all-positive
    MREgger.res4 <- MREgger(gamma_hat, Gamma_hat, sigma_X, sigma_Y, weighting=1, recode=F, remove_weak=F)[2,] #First-order weights, estimated SNP-X coefficients, no recoding
    MREgger.res5 <- MREgger(gamma_hat, Gamma_hat, sigma_X, sigma_Y, weighting=1, recode=T, remove_weak=T)[2,] #First-order weights, estimated SNP-X coefficients, recode to all-positive
    MREgger.res6 <- MREgger(gamma_hat, Gamma_hat, sigma_X, sigma_Y, weighting=0, recode=T, remove_weak=F)[2,] #Unweighted, estimated SNP-X coefficients, recode to all-positive
    
    MREgger.res  <- rbind(MREgger.res1,  MREgger.res2,  MREgger.res3, MREgger.res4,  MREgger.res5,  MREgger.res6)
    
    #Heterogeneity test
    Het.res1 <- Het_test(gamma_hat, Gamma_hat, sigma_X, sigma_Y, weighting=1, recode=T, remove_weak=F)[c(4,8)] #First-order weights, estimated SNP-X coefficients, recode to all-positive
    Het.res2 <- Het_test(gamma_hat, Gamma_hat, sigma_X, sigma_Y, weighting=2, recode=T, remove_weak=F)[c(4,8)] #Modified weights, estimated SNP-X coefficients, recode to all-positive
    Het.res3 <- Het_test(gamma,     Gamma_hat, sigma_X, sigma_Y, weighting=1, recode=T, remove_weak=F)[c(4,8)] #First-order weights, true SNP-X coefficients, recode to all-positive
    Het.res4 <- Het_test(gamma_hat, Gamma_hat, sigma_X, sigma_Y, weighting=1, recode=F, remove_weak=F)[c(4,8)] #First-order weights, true SNP-X coefficients, no recoding
    Het.res5 <- Het_test(gamma_hat, Gamma_hat, sigma_X, sigma_Y, weighting=1, recode=T, remove_weak=T)[c(4,8)] #First-order weights, estimated SNP-X coefficients, recode to all-positive

    Het.res  <- c(Het.res1, Het.res2, Het.res3, Het.res4, Het.res5)
    
    #Heteroscedasticity test
    Var.res1 <- Var_test(gamma_hat, Gamma_hat, sigma_X, sigma_Y, weighting=1, recode=T, remove_weak=F) #First-order weights, estimated SNP-X coefficients, recode to all-positive
    Var.res2 <- Var_test(gamma_hat, Gamma_hat, sigma_X, sigma_Y, weighting=2, recode=T, remove_weak=F) #Modified weights, estimated SNP-X coefficients, recode to all-positive
    Var.res3 <- Var_test(gamma,     Gamma_hat, sigma_X, sigma_Y, weighting=1, recode=T, remove_weak=F) #First-order weights, true SNP-X coefficients, recode to all-positive
    Var.res4 <- Var_test(gamma_hat, Gamma_hat, sigma_X, sigma_Y, weighting=1, recode=F, remove_weak=F) #First-order weights, estimated SNP-X coefficients, no recoding
    Var.res5 <- Var_test(gamma_hat, Gamma_hat, sigma_X, sigma_Y, weighting=1, recode=T, remove_weak=T) #First-order weights, estimated SNP-X coefficients, no recoding
    Var.res6 <- Var_test(gamma_hat, Gamma_hat, sigma_X, sigma_Y, weighting=0, recode=T, remove_weak=F) #Unweighted, estimated SNP-X coefficients, recode to all-positive
    
    Var.res  <- c(Var.res1, Var.res2, Var.res3, Var.res4, Var.res5, Var.res6)
    
    c(as.numeric(apply(MREgger.res, 1, FUN=get.res, beta)),
      Het.res<0.05,
      Var.res<0.05)
  }  
  sim.res <- rbind(sim.res, c(apply(sim.res.a, 2, mean)))
}

#Format the results
sim.res <- data.frame(scenarios, sim.res)

colnames_MREgger <- paste(paste('MREgger', rep(1:6,each=4), sep=''), rep(c('bias','se','coverage','rr'),6), sep='_')
colnames_HetTest <- paste(rep(c('Q','Qline'),5), rep(1:5,each=2), sep='_')
colnames_VarTest <- paste(rep(c('Linear','Quadratic'),6), rep(1:6,each=2), sep='_')

colnames(res)[8:ncol(res)] <- c(colnames_MREgger, colnames_HetTest, colnames_VarTest)

write.table(sim.res, 'Simulation_results.txt', row.names=F, sep='\t')