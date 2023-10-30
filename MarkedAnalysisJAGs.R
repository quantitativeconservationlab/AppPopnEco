######################################################################
#######################################################################
####   This script was created by Dr. Jen Cruz to estimate     ########
# abundance of ground squirrels at the NCA from capture-recapture     #
# surveys conducted at 20 siteXyear during 2021-2023.                 #
# Detection j predictors: 1-effort: number of hours the traps were open #
# each day, 2-mean tempC_st up to the hour the traps were closed, 3-mean #
# wind up to the hour the traps were closed each day.                 #
#                                                                     #
# Predictors for abundance: shrub, annual, perennial, herbaceous      #

#Data cleaned in MarkedPrep.R                                         #
#######################################################################
##### Set up your workspace and load relevant packages -----------
# load packages:
library( tidyverse ) 
#library( sf )
library( jagsUI )
## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are. 
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# load workspace with clean data
load( "MarkPrepWorkspace.RData")
###########################################################################
####################### define MCMC settings ##############################

ni <- 20000; nt <- 20; nb <- 100000; nc <- 5 #iterations, thinning, burnin, chains

##### end of MCMC parameters definition ##############
###########################################################################
###### m1 single season abundance model with NO data augmentation  #
# y[i]|N[i]~Multinom(N[i], pi ) #
# from chpt 7 pg342 of applied hierarchical model book #
# N is no longer in the model but we can still obtain predictions of it #
# for every site by N[i]~Poisson(lambda[i]). This is not conditional on #
#observed individuals but are predictions at sites as if they were not sampled #
# we also add extra residual overdispersion using Poison/loglinear version#
# by adding a random effect for site
# ecological predictors: int + eps[i] + habitat[i] 
# detection predictors: wind + temp + effort 
############################################################################
############## Specify model in bugs language:  #####################
sink( "m1.txt" )
cat( "
     model{
     
      #priors
      #for detection model: 
      #define intercept as mean probs:
      int.det <- log( mean.det / ( 1 - mean.det ) )
      mean.det ~ dbeta( 4, 4 )
      
      #priors for detection coefficients:
      for( a in 1:A ){
        #define as a slightly informative prior
        alpha[a] ~ dnorm( 0, 0.2 ) T(-7, 7 )
      }
      
      #define intercept for abundance just with normal
      int.lam ~ dnorm( 0, 0.01 )

      #priors for abundance coefficients:
      for( b in 1:B ){
        #define as a slightly informative prior
      beta[b] ~ dnorm( 0, 0.2 ) T(-7, 7 )
       }
      
    # Convert to Poisson/lognormal to account for overdispersion 
      # between sites by adding a random intercept for site
      # #random intercept for site #did not converge 
      prec.i <- 1 / ( sigma.i * sigma.i )
      sigma.i ~ dt( 0, 2.5, 7 ) T(0, )
      
    #random site intercept
    # makes the model a Poisson/log Normal accommodating 
    # for overdispersion
      for ( i in 1:S ){  #loop over sites
        eps.i[ i ] ~ dnorm( 0, prec.i )
      } #S

    # log-linear ecological model of abundance
    for( i in 1:I ){ #sitesXyear combo

      #log-linear model for abundance model 
      log( lambda[ i ] ) <- int.lam + 
                        eps.i[ site[i] ] + #random site effect
                        #we matrix multiply for all habitat predictors:
                        inprod( beta, XIK[ i, ])
        
          #poisson parameter = multinomial cell probabilities of 
          # expected abundance for each survey day
          pi[i,1] <- (1-p[i,1]) * (1 - p[i,2]) * p[i,3] * lambda[ i ] #001
          pi[i,2] <- (1-p[i,1]) * p[i,2]  * (1-p[i,3]) * lambda[ i ]#010
          pi[i,3] <- (1-p[i,1]) * p[i,2] * p[i,3] * lambda[ i ] #011
          pi[i,4] <- p[i,1] * (1-p[i,2]) * (1-p[i,3]) * lambda[ i ] #100
          pi[i,5] <- p[i,1] * (1-p[i,2]) * p[i,3] * lambda[ i ] #101
          pi[i,6] <-  p[i,1] * p[i,2] * (1-p[i,3]) * lambda[ i ] #110
          pi[i,7] <- p[i,1] * p[i,2] *p[i,3] * lambda[ i ] #111
      
      #loop over each column      
      for( c in 1:CH ){
      
        #relate counts of capture histories to detectionXlambda
        y[ i, c ] ~ dpois( pi[i,c] )
        
        #predict capture histories
        y_hat[ i, c ] ~ dpois( pi[i,c] )
        #calculate model residuals
        e1[ i,c ] <- pi[ i,c ] * N[i]
        resid1[i,c] <- pow( pow( y[ i, c ], 0.5 ) - 
                          pow( e1[i,c], 0.5), 2 )
        resid1_hat[i,c] <- pow( pow( y_hat[ i, c ], 0.5 ) - 
                          pow( e1[i,c], 0.5), 2 )
        
       }
        
      for( j in 1:J ){ #loop over survey days
      
        logit( p[ i, j ] ) <- int.det +
                  #fixed effects
                  alpha[1] * wind[ i, j ] +
                  alpha[2] * temp[ i, j ] +
                  alpha[3] * effort[ i, j ]
        
      } #close J
      
      #generate predictions of N[i]
      N[ i ] ~ dpois( lambda[i] )
      
      }#close i

      #fit statistic for observed data y, 
      #think of it as fit of encounter model
      #how well does the model for y|n fit?
      fit1 <- sum( resid1[,] )
      fit1_hat <- sum( resid1_hat[,] )

     } #model close
     
     ", fill = TRUE )

sink()

################ end of model specification  #####################################     
modelname <- "m1.txt"
#parameters monitored
params <- c(  'int.det' #intercept for detection
              , 'int.lam' #intercept for lambda
              , 'alpha' #detection coefficients
              , 'beta' #abundance coefficients
              , 'eps.i' #random effect for site
              , 'sigma.i' #error for site random intercept
              , 'N' #estimates of abundance by siteXyear
              , 'p' #detection
              , 'y_hat' #predicted capture histories
              , "fit1", "fit1_hat" #fit statistics
)

head( ik_sc)
# We select the abundance predictors we want
XIK <- ik_sc[,c("shrub", "perennial", "annual")]
#how many ecological predictors
B <- dim(XIK)[2]
#how many detection predictors
A <- 3
#define initial parameter values
inits <- function(){ list( beta = rnorm( B ),
                           alpha = rnorm( A ),
                           int.lam = runif(1), int.p = runif(1) ) }

#define data that will go in the model
#note that we do not have to remove sites with all zero captures here:
str( win.data <- list( y = as.matrix(y_ik),
                       #number of sitesXyear, surveys, det predictors, and abund preds
                       I = I, J = J, A = A, B = B, CH = CH,
                       #site id for each unique site
                       site = as.numeric(as.factor(ik_df$SiteID)),
                       #number of sites
                       S = length( unique(ik_df$SiteID)),
                       #site level habitat predictors
                       XIK = XIK,
                       #observation predictors only for sites with captures
                       wind = ij_wide[,widx],
                       temp = ij_wide[,tidx],
                       effort = ij_wide[,eidx]
) )

#call JAGS and summarize posteriors:
m1 <- autojags( win.data, inits = inits, params, modelname, #
                  n.chains = nc, n.thin = nt, n.burnin = nb,
                  iter.increment = ni, max.iter = 1000000, 
                  Rhat.limit = 1.05,
                  save.all.iter = FALSE, parallel = TRUE ) 
#m1 <- update( m1, parameters.to.save= params,
#                    n.iter = 500000, n.thin = nt)
#
plot(m1)
summary(m1)

###### end m1 ########
###########################################################################
###### m2 single season abundance model with NO data augmentation  #
# from chpt 7 pg334 of applied hierarchical model book, option 3,   #
# conditional multinomial. This method breaks multinomial into 2: #
# 1. part conditioned on known sample size (individuals encountered), #
# 2. binomial with unknown sample size                                # 
#  N[i]~Poisson(lambda[i]); n[i] ~ Binomial( N[i], 1 - pi[0] );      #
# y[i]|n[i] ~ Multi( n[i], pi[i] )  #
# ecological predictors: int + eps[i] #random intercept for site    #
# detection predictors: wind + temp + effort                     #
############################################################################
############## Specify model in bugs language:  #####################
sink( "m2.txt" )
cat( "
     model{
     
      #priors
      #for detection model: 
      #define intercept as mean probs:
      int.det <- log( mean.det / ( 1 - mean.det ) )
      mean.det ~ dbeta( 4, 4 )
      
      #priors for detection coefficients:
      for( a in 1:A ){
        #define as a slightly informative prior
        alpha[a] ~ dnorm( 0, 0.2 ) T(-7, 7 )
      }
      
      #define intercept for abundance just with normal
      int.lam ~ dnorm( 0, 0.01 )

      #priors for abundance coefficients:
      for( b in 1:B ){
        #define as a slightly informative prior
        beta[b] ~ dnorm( 0, 0.2 ) T(-7, 7 )
       }
      
    # Convert to Poisson/lognormal to account for overdispersion 
      # between sites by adding a random intercept for site
      # #random intercept for site #did not converge 
      prec.i <- 1 / ( sigma.i * sigma.i )
      sigma.i ~ dt( 0, 2.5, 7 ) T(0, )
      
    #random site intercept
    # makes the model a Poisson/log Normal accommodating 
    # for overdispersion
      for ( i in 1:S ){  #loop over sites
        eps.i[ i ] ~ dnorm(0, prec.i )
      } #M


    # log-linear ecological model of abundance
    for( i in 1:I ){#sitesXyear combo

      #log-linear model for abundance model 
      log( lambda[ i ] ) <- int.lam + 
                        eps.i[ site[i] ] + #random site effect
                        #effect of herbaceous cover
                        inprod( beta, XIK[ i, ] )
        
          #poisson parameter = multinomial cell probabilities of 
          # expected abundance for each survey day
          pi[i,1] <- (1-p[i,1]) * (1 - p[i,2]) * p[i,3] #001
          pi[i,2] <- (1-p[i,1]) * p[i,2]  * (1-p[i,3]) #010
          pi[i,3] <- (1-p[i,1]) * p[i,2] * p[i,3]  #011
          pi[i,4] <- p[i,1] * (1-p[i,2]) * (1-p[i,3])#100
          pi[i,5] <- p[i,1] * (1-p[i,2]) * p[i,3]  #101
          pi[i,6] <-  p[i,1] * p[i,2] * (1-p[i,3])  #110
          pi[i,7] <- p[i,1] * p[i,2] *p[i,3] #111
        #probability of 000 capture history
        pi0[i] <- 1 - (pi[i,1] + pi[i,2] +pi[i,3] +pi[i,4] +pi[i,5]+
                    pi[i,6] + pi[i,7] )
        pcap[i] <- 1 - pi0[i]
          
      #loop over each column      
      for( c in 1:CH ){
        pic[ i, c ] <- pi[ i, c ] / pcap[ i ]
        
        #estimated y based on encounter model
        e1[i,c] <- pic[ i, c ] * n[i]
        resid1[i,c] <- pow( pow( y[ i, c ], 0.5 ) - pow( e1[i,c], 0.5), 2 )
        resid1_hat[i,c] <- pow( pow( y_hat[ i, c ], 0.5 ) - 
                          pow( e1[i,c], 0.5), 2 )
       }
        
      #observations based on conditional cell probabilities  
      y[ i, 1:CH ] ~ dmulti( pic[ i, 1:CH ], n[ i ] )  
      #estimated capture histories
      y_hat[ i, 1:CH ] ~ dmulti( pic[ i, 1:CH ], n[ i ] )  
      
      #model for observed sample size
      n[ i ] ~ dbin( pcap[ i ], N[ i ] )
      #estimated number of indiv. observed
      n_hat[ i ] ~ dbin( pcap[ i ], N[ i ] )
      
      #expected n based on model
      e2[ i ] <- pcap[ i ] * lambda[ i ]
      resid2[ i ] <-  pow( pow( n[ i ], 0.5 ) - pow( e2[ i ], 0.5), 2 )
      resid2_hat[ i ] <- pow( pow( n_hat[ i ], 0.5 ) - pow( e2[ i ], 0.5), 2 )
      
      #process model
      N[ i ] ~ dpois( lambda[i] )

      for( j in 1:J ){ #loop over survey days
      
        logit( p[ i, j ] ) <- int.det +
                  #fixed effects
                  alpha[1] * wind[ i, j ] +
                  alpha[2] * temp[ i, j ] +
                  alpha[3] * effort[ i, j ]
        
      } #close J
      
      }#close n
    
    #Bayesian pvalue based on freeman-tukey statistic where 
    # FIT = sum( y^0.5 - e^0.5)^2 
    
    #fit statistic for observed data y, think of it as fit of encounter model
    #how well does the model for y|n fit?
    fit1 <- sum( resid1[,] )
    fit1_hat <- sum( resid1_hat[,] )
    
    #fit of model for n. Think of it as a test for the abundance part of the model #
    #since having variation in N shows up in the expected value of n.
    fit2 <- sum( resid2[] )
    fit2_hat <- sum( resid2_hat[] )
    
     } #model close
     
     ", fill = TRUE )

sink()

################ end of model specification  #####################################     
modelname <- "m2.txt"
#parameters monitored
params <- c(  'int.det' #intercept for detection
              , 'int.lam' #intercept for lamda
              , 'alpha' #detection coefficients
              , 'beta' #abundance coefficients
              , 'eps.i' #random effect for site
              , 'sigma.i' #error for site random intercept
              , 'N' #estimates of abundance by siteXyear
              , 'p' #detection
              , 'y_hat' #predicted y
              , "fit1", "fit1_hat", "fit2", "fit2_hat" #fit statistics
)

#define initial parameter values
inits <- function(){ list( beta = rnorm( B ),
                           alpha = rnorm( A ),
                           int.lam = runif(1), int.p = runif(1),
                           N = n_ik + 2 ) }

#define data that will go in the model
str( win.data <- list( y = as.matrix(y_ik), 
                       #number of trapped individuals at each site
                       n = n_ik,
                       #number of sitesXyear, surveys, det predictors, and abund preds
                       I = I, J = J, A = A, B = B, CH = CH,
                       #site id for each unique site
                       site = as.numeric(as.factor(ik_df$SiteID)),
                       #number of sites
                       S = length( unique(ik_df$SiteID)),
                       #site level habitat predictors
                       XIK = XIK,
                       #observation predictors only for sites with captures
                       wind = ij_wide[,widx],
                       temp = ij_wide[,tidx],
                       effort = ij_wide[,eidx]
) )

#call JAGS and summarize posteriors:
m2 <- autojags( win.data, inits = inits, params, modelname, #
                  n.chains = nc, n.thin = nt, n.burnin = nb,
                  iter.increment = ni, max.iter = 1000000, 
                  Rhat.limit = 1.05,
                  save.all.iter = FALSE, parallel = TRUE ) 
#m2 <- update( m2, parameters.to.save= params,
#                    n.iter = 500000, n.thin = nt)
#
plot(m2)
summary(m2)

###### end m2 ########


##### save relevant stuff ##################################
save.image( "MarkedResultsJAGs.RData")

################### end of script #######################################
