rm(list=ls())

dirs <- c('C:/Users/Fernando Hartwig/Dropbox/Methodological papers/InSIDE_test/PLOS Genetics/',
          'C:/Users/ICEH_DESK12/Dropbox/Methodological papers/InSIDE_test/PLOS Genetics/')
dir <- dirs[dir.exists(dirs)]
setwd(dir)

#Install (if needed) and load the lmtest package
if(!requireNamespace('lmtest', quietly=TRUE)) { install.packages('lmtest') }; library(lmtest)

##########################
# User-defined functions #
##########################

#Function for nicely-formated rouding
rounding <- function(input, digits=1) {
  
  output <- NULL
  
  for(x in input) {
    if(!is.na(x)) {
      
      #Round the numbers and conver to character
      x <- as.character(round(x, digits=digits))
      
      ind_e <- grepl('e', x, fixed=T)
      if(ind_e) {
        x.split <- strsplit(x, split='e', fixed=T)[[1]]
        x       <- as.character(round(as.numeric(x.split[1]), digits=digits))
      }
      
      #Identify negative numbers, stored this info and remove the minus sign
      neg.x <- substr(x, 1, 1)=='-'
      x[neg.x] <- substr(x[neg.x], 2, nchar(x[neg.x])) 
      
      #Identify numbers without decimals
      decimal.split <- strsplit(x, split='.', fixed=T)
      length.after.decimal <- numeric(length(decimal.split))
      
      for(a in 1:length(decimal.split)) {
        cur.split <- decimal.split[[a]]
        if(length(cur.split)==2) {
          length.after.decimal[a] <- nchar(cur.split[2]) 
        }
      }
      
      #Add zeros
      for(a in 1:length(x)) {
        cur.x <- x[a]
        
        if(length.after.decimal[a]!=digits) {
          if(length.after.decimal[a]==0) {
            x[a] <- paste(cur.x, '.',
                          paste(rep('0', digits-length.after.decimal[a]), collapse=''),
                          sep='')
          } else {
            x[a] <- paste(cur.x,
                          paste(rep('0', digits-length.after.decimal[a]), collapse=''),
                          sep='')
          }
        }
      }
      
      #Re-insert negative sign
      x[neg.x] <- paste('-',x[neg.x], sep='')
      
      if(ind_e) {
        x <- paste(x, 'e', x.split[2], sep='')
      }
    }
    output <- c(output, x)
  }
  return(output)
}

#Inverse variance weighted estimator
IVW <- function(BetaXG,        #BetaXG:   SNP-exposure coefficient, typically a linear (if the exposure is continuous) or logistic (if the exposure is binary) regression coefficient.
                BetaYG,        #BetaYG:   SNP-outcome coefficient, typically a linear (if the outcome is continuous) or logistic (if the outcome is binary) regression coefficient.
                seBetaXG,      #seBetaXG: standard error of BetaXG.
                seBetaYG,      #seBetaYG: standard error of BetaYG.
                alpha=0.05) {  #alpha:    1-alpha is the confidence level of confidence intervals.
  
  #IVW using 1st order weights
  IVWfit  <- summary(lm(BetaYG ~ -1+BetaXG, weights=seBetaYG^-2))
  
  #Inference with correct standard errors  
  DF      <- length(BetaYG)-1
  IVWBeta <- IVWfit$coef[1,1]
  SE      <- IVWfit$coef[1,2]/min(IVWfit$sigma, 1)
  IVW_p   <- 2*(1-pt(abs(IVWBeta/SE),DF))
  IVW_CI  <- IVWBeta + c(-1,1)*qt(df=DF, 1-alpha/2)*SE
  
  #IVWResults = (point estimate, corrected standard error, 
  #95% Confidence interval, t-statistic, p-value)  
  IVWResults <- c(IVWBeta,SE,IVW_CI,IVWBeta/SE,IVW_p)
  names(IVWResults) <- c('Estimate', 'SE', 'CI_low', 'CI_upp', 'T', 'P-value')  
  
  return(IVWResults)
}

#MR-Egger regression
MREgger <- function(BetaXG, BetaYG, seBetaXG, seBetaYG, alpha=0.05) {

  #Pre-processing steps to ensure all instrument-exposure estimates are positive  
  BYG <- BetaYG*sign(BetaXG)
  BXG <- abs(BetaXG)  

  #MR-Egger model with 1st order weights
  MREggerFit <- summary(lm(BYG~BXG, weights=seBetaYG^-2))
  
  #Inference with correct standard errors
  MREggerBeta0 <- MREggerFit$coef[1,1]
  MREggerBeta1 <- MREggerFit$coef[2,1]
  SE0          <- MREggerFit$coef[1,2]/min(MREggerFit$sigma, 1)
  SE1          <- MREggerFit$coef[2,2]/min(MREggerFit$sigma, 1)
  DF           <- length(BetaYG)-2
  MRBeta0_p    <- 2*(1-pt(abs(MREggerBeta0/SE0),DF))
  MRBeta1_p    <- 2*(1-pt(abs(MREggerBeta1/SE1),DF))
  MRBeta0_CI   <- MREggerBeta0 + c(-1,1)*qt(df=DF, 1-alpha/2)*SE0
  MRBeta1_CI   <- MREggerBeta1 + c(-1,1)*qt(df=DF, 1-alpha/2)*SE1
  
  # MREggerResults = (point estimate, corrected standard error, 
  # 95% Confidence interval, t-statistic, p-value) for
  # intercept (row 1) and slope (row 2). 
  MREggerResults <- matrix(nrow = 2,ncol = 6)
  MREggerResults[1,] <- c(MREggerBeta0,SE0,MRBeta0_CI,MREggerBeta0/SE0,MRBeta0_p)  
  MREggerResults[2,] <- c(MREggerBeta1,SE1,MRBeta1_CI,MREggerBeta1/SE1,MRBeta1_p)
  colnames(MREggerResults) <- c('Estimate', 'SE', 'CI_low', 'CI_upp', 'T', 'P-value')
  rownames(MREggerResults) <- c('Intercept', 'Beta')  
  
  return(MREggerResults)
}

#Median estimator
MedianMR <- function(BetaXG, BetaYG, seBetaXG, seBetaYG, alpha=0.05,
                     n_boot=1e4) { #n_boot: number of bootstrap replicates (used to calculate standard errors)
  
  #Define a function to obtain the WM estimate
  weighted.median <- function(betaIV.in, weights.in) {
    betaIV.order <- betaIV.in[order(betaIV.in)]
    weights.order <- weights.in[order(betaIV.in)]
    weights.sum <- cumsum(weights.order)-0.5*weights.order
    weights.sum <- weights.sum/sum(weights.order)
    below <- max(which(weights.sum<0.5))
    weighted.est <- betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
      (0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
    return(weighted.est)
  }
  
  #Define a function to calculate standard errors through boostraping
  median.boot <- function(betaXG.in, betaYG.in, sebetaXG.in, sebetaYG.in, weights.in) {
    med <- NULL
    for(i in 1:n_boot) {
      betaXG.boot <- rnorm(length(betaXG.in), mean=betaXG.in, sd=sebetaXG.in)
      betaYG.boot <- rnorm(length(betaYG.in), mean=betaYG.in, sd=sebetaYG.in)
      betaIV.boot <- betaYG.boot/betaXG.boot
      med[i] <- weighted.median(betaIV.boot, weights.in)
    }
    return(sd(med))
  }
  
  betaIV <- BetaYG/BetaXG #ratio estimates
  weights <- (seBetaYG/BetaXG)^-2 #1st order weights
  
  #Calculate point estimates and associated standard errors
  #SM: simple median. WM: weighted median.
  betaSM <- median(betaIV)
  betaWM <- weighted.median(betaIV, weights)
  
  sebetaSM <- median.boot(BetaXG, BetaYG, seBetaXG, seBetaYG, rep(1, length(BetaXG))/length(BetaXG))  
  sebetaWM <- median.boot(BetaXG, BetaYG, seBetaXG, seBetaYG, weights)
  
  #Calculate (1-alpha)*100% confidence intervals (CIs) using the normal approximation
  betaSM_CI <- betaSM + c(-1,1)*qnorm(1-alpha/2)*sebetaSM
  betaWM_CI <- betaWM + c(-1,1)*qnorm(1-alpha/2)*sebetaWM
  
  #Calculate P-value
  pSM <- pt(abs(betaSM/sebetaSM), df=length(BetaXG)-1, lower.tail=F)*2
  pWM <- pt(abs(betaWM/sebetaWM), df=length(BetaXG)-1, lower.tail=F)*2
  
  ResultsSM <- c(betaSM, sebetaSM, betaSM_CI, betaSM/sebetaSM, pSM)
  ResultsWM <- c(betaWM, sebetaWM, betaWM_CI, betaWM/sebetaWM, pWM)  
  
  Results <- rbind(ResultsSM, ResultsWM)
  
  colnames(Results) <- c('Estimate', 'SE', 'CI_low', 'CI_upp', 'T', 'P-value')
  rownames(Results) <- c('Simple median', 'Weighted median')
  
  return(Results)
}

#Heterogeneity test
Het_test <- function(BetaXG, BetaYG, seBetaXG, seBetaYG) {
  
  #Cocharn's Q
  BetaIV      <- BetaYG/BetaXG                   #Ratio estimates
  weights1    <- (seBetaYG/abs(BetaXG))^-2       #1st order weights
  IVW_beta_w1 <- weighted.mean(BetaIV, weights1) #IVW point estimate using 1st order weights
  
  weights2 <- 1/(seBetaYG^2/BetaXG^2 + (IVW_beta_w1^2)*seBetaXG^2/BetaXG^2) #Modified 2nd order weights
  BIVw2    <- BetaIV*sqrt(weights2)
  sW2      <- sqrt(weights2)
  IVWfit   <- summary(lm(BIVw2 ~ -1+sW2)) #IVW model with modified 2nd order weights
  
  DF         <- length(BetaIV)-1
  phi_IVWfit <- IVWfit$sigma^2
  Q_stat     <- DF*phi_IVWfit                  #Q statistic
  Q_I2       <- max(0, (Q_stat-DF)/Q_stat)*100 #I2
  Q_p        <- 1-pchisq(Q_stat, DF)           #P-value
  
  #Rucker's Q'
  BYG <- BetaYG*sign(BetaXG) # Pre-processing steps to ensure all  
  BXG <- abs(BetaXG)         # gene-exposure estimates are positive  
  
  weights1        <- 1/seBetaYG^2 #1st order weights
  MREggerBeta1_w1 <- coef(lm(BYG ~ BXG, weights=weights1))[2] #MR-Egger point estimate with 1st order weights
  
  weights2   <- 1/(seBetaYG^2 + (MREggerBeta1_w1^2)*seBetaXG^2) #Modified 2nd order weights
  MREggerFit <- summary(lm(BYG ~ BXG, weights=weights2)) #MR-Egger model with modified 2nd order weights
  
  DF             <- length(BXG)-2
  phi_MREggerfit <- MREggerFit$sigma^2
  Qline_stat     <- DF*phi_MREggerfit                      #Q' statistic
  Qline_I2       <- max(0, (Qline_stat-DF)/Qline_stat)*100 #I2
  Qline_p        <- 1-pchisq(Qline_stat, DF)               #P-value
  
  return(c(Q_stat,     DF+1, Q_I2,     Q_p,
           Qline_stat, DF,   Qline_I2, Qline_p))
}

#Heteroscedasticity test
Var_test <- function(BetaXG, BetaYG, seBetaXG, seBetaYG) {
  
  #Calculate weights
  weights <- 1/seBetaYG^2 #1st order weights
  
  #Pre-processing steps to ensure all instrument-exposure estimates are positive
  BYG <- BetaYG*sign(BetaXG)
  BXG <- abs(BetaXG)  
  
  weights <- length(weights)*weights/sum(weights) #Normalize weights so they sum to the number of genetic variants
  
  #Heteroscedasticity test
  fit <- lm(BYG ~ BXG, weights=weights) #MR-Egger model
  
  pval <- bptest(fit, varformula = ~ BXG)$p.value

  return(pval)
}

###############
# MR analyses #
###############

#Load data from supplementary material. The code below assumes the data is in a dataframe named "data" containing four columns: BetaXG, seBetaXG, BetaYG, seBetaYG
res <- NULL

#Two sets of analyses:
#a=1: uses all 27 variants
#a=2: excludes the single clear outlier
for(a in 1:2) {
  
  if(a==1) {
    ind <- rep(T, nrow(data))
    
  } else {
    ind[21] <- F
  }
  
  data.a <- data[ind,]
  
  #MR analyses
  IVW.res     <- IVW(BetaXG=data.a$BetaXG, BetaYG=data.a$BetaYG, seBetaXG=data.a$seBetaXG, seBetaYG=data.a$seBetaYG)[c(1,3,4)]
  Median.res  <- MedianMR(BetaXG=data.a$BetaXG, BetaYG=data.a$BetaYG, seBetaXG=data.a$seBetaXG, seBetaYG=data.a$seBetaYG)[,c(1,3,4)]
  MREgger.res <- MREgger(BetaXG=data.a$BetaXG, BetaYG=data.a$BetaYG, seBetaXG=data.a$seBetaXG, seBetaYG=data.a$seBetaYG)[-1,c(1,3:4)]
  res.MR      <- rbind(IVW.res, Median.res, MREgger.res)
  
  res.a <- paste(rounding(res.MR[,1],2), ' (',
                 rounding(res.MR[,2],2), ' to ',
                 rounding(res.MR[,3],2), ')', sep='')
  
  #Heterogeneity tests
  Het.res <- Het_test(BetaXG=data.a$BetaXG, BetaYG=data.a$BetaYG, seBetaXG=data.a$seBetaXG, seBetaYG=data.a$seBetaYG)
  
  cochrans.q <- Het.res[c(1,2,4)]
  p          <- cochrans.q[3]
  if(p<0.0001) { p <- 'P<0.0001'} else { p <- paste('P=', rounding(p,4), sep='')  }
  cochrans.q <- paste(rounding(cochrans.q[1],2), ' (df=', cochrans.q[2], '); ', p, sep='')
  
  ruckers.qline <- Het.res[c(5,6,8)]
  p          <- ruckers.qline[3]
  if(p<0.0001) { p <- 'P<0.0001'} else { p <- paste('P=', rounding(p,4), sep='')  }
  ruckers.qline <- paste(rounding(ruckers.qline[1],2), ' (df=', ruckers.qline[2], '); ', p, sep='')
  
  res.a <- c(res.a, cochrans.q, ruckers.qline)
  
  #Heteroscedasticity tests
  Var.res <- Var_test(BetaXG=data.a$BetaXG, BetaYG=data.a$BetaYG, seBetaXG=data.a$seBetaXG, seBetaYG=data.a$seBetaYG)
  
  if(Var.res<0.0001) { Var.res <- 'P<0.0001'} else { Var.res <- paste('P=', rounding(Var.res,4), sep='')  }

  res.a <- c(res.a, Var.res)
  
  #Correlation between the absolute value of MR-Egger residuals and instrument strength, as a further metric to quantify heteroscedasticity
  
  #Pre-processing steps to ensure all instrument-exposure estimates are positive
  BYG <- data.a$BetaYG*sign(data.a$BetaXG)
  BXG <- abs(data.a$BetaXG)
  
  #Calculate weights
  weights1        <- 1/data.a$seBetaYG^2                                           #1st order weights
  MREggerBeta1_w1 <- coef(lm(BYG ~ BXG, weights=weights1))[2]                      #MR-Egger point estimate with 1st order weights
  weights2        <- 1/(data.a$seBetaYG^2 + (MREggerBeta1_w1^2)*data.a$seBetaXG^2) #Modified 2nd order weights
  fit             <- lm(BYG ~ BXG, weights=weights2) #MR-Egger model with modified 2nd order weights
  cor.res         <- rounding(cor(abs(fit$residuals),BXG),2)
  
  res.a <- c(res.a, cor.res)
  
  res <- cbind(res, res.a)
}

#Print results on screen
for(a in 1:nrow(res)) {
  cat(res[a,], sep='\t'); cat('', sep='\n')
}

#Check if observed heteroscedasticity is attributable to variability in seBetaYG according to BXG
n.sim <- 1e4

p.het.sim <- numeric(n.sim)

for(a in 1:n.sim) {
  
  BetaYG.sim   <- rnorm(nrow(data), 0, sample(data$seBetaYG))
  p.het.sim[a] <- Var_test(BetaXG=data$BetaXG, BetaYG=BetaYG.sim, seBetaXG=data$seBetaXG, seBetaYG=data$seBetaYG)
  
}

mean(p.het.sim)
mean(p.het.sim<0.05)