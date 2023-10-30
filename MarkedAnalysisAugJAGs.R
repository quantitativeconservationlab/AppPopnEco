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
#########################################################################
###### ma1 single season abundance model using trapping with augmentation#
# code from Chpt 7 section 7.8.4 from applied hierarchical modeling book#
#  Mht model including random indiv intercepts for individuals          #
# ecological predictors: lambda[i] ~ int + hab[i] + eps[i]              #
# detection predictors: p[n,j,i] ~ int + wind[i,j] + temp[i,j] +              #
#  effort[i,j] + eta[n] #random individual intercept could be replaced  #
# with individual covariates                                            #
#########################################################################
############## Specify model in bugs language:  #####################
sink( "ma1.txt" )
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
        alpha[ a ] ~ dnorm( 0, 0.2 ) T(-7, 7 )
      }
      
      #define intercept for abundance just with normal
      int.lam ~ dnorm( 0, 0.01 )

      # #priors for abundance coefficients:
      #for( b in 1:B ){
        #define as a slightly informative prior
       beta ~ dnorm( 0, 0.2 ) T(-7, 7 )
       #}
      
      #prior for individual random effects 
       prec.n <- 1 / ( sigma.n * sigma.n )
       sigma.n ~ dt( 0, 2.5, 7 ) T(0, )
      
      # #random intercept for site #did not converge 
      prec.i <- 1 / ( sigma.i * sigma.i )
      sigma.i ~ dt( 0, 2.5, 7 ) T(0, )
      
      #estimated prob that animal was present
      #M is number of augmented individuals, 
      #lambda is summed across all sites sum(lambda[1:I])
      psi <- sum( lambda[]) / M 

    # log-linear ecological model of abundance
    for( i in 1:I ){#sites
      # #random intercept for site
      # which makes the model a Poisson/log Normal accommodating 
      # for overdispersion (alternatives include negative binomial or ZIP)
      eps[i] ~ dnorm(0, prec.i )
      
      #probs of belonging to a site
      probs[i] <- lambda[i] / sum( lambda[] )
      
      log( lambda[ i ] ) <- int.lam +
                          #random site effect    
                          eps[i] +
                        inprod( beta, XIK[ i ]  )
                        
    } #i
    
      #model for individual encounter histories:
      for( n in 1:M ){ #loop over all individuals
      
        #site membership
        site[ n ] ~ dcat( probs[] )
        
        #was animal present and undetected or not present
        z[ n ] ~ dbern( psi )
        
        #individual random effect
        eta[ n ] ~ dnorm( 0, prec.n ) 
        
        for( j in 1:J ){ #loop over survey occassions
        
          logit( p[n,j] ) <-  int.det + eta[ n ] + 
                  #fixed effects
                  alpha[1] * wind[ site[n], j ] +
                  alpha[2] * temp[ site[n], j ] +
                  alpha[3] * effort[ site[n], j ]
          
          #relate observed counts to abundance and detection
          pz[ n, j ] <- p[ n, j ] * z[ n ]
          y[ n, j ] ~ dbern( pz[n,j] )  
        
      } #close J
      }#close n

      #derived parameters
      for( n in 1:M ){
        site.out[ n ] <- site[ n ] * z[ n ]

        for( i in 1:I ){#sites
          #estimated abundance
          N.site[i,n] <- step( 0.01 * ( i - site.out[n] ) - 0.02 *
                ( i - site.out[n] ) * ( i - site.out[n] ) + 0.001 )
      }}

    N.tot <- sum( z[1:M] )
    
     } #model close
     
     ", fill = TRUE )

sink()
################ end of model specification  #####################################     
modelname <- "ma1.txt"
#parameters monitored
params <- c(  'int.det' #intercept for detection
              , 'int.lam' #intercept for lamda
              , 'alpha' #detection coefficients
              , 'sigma.n' #error for random intercept
              , 'beta' #abundance coefficients
              , 'sigma.i', 'sigma.n' #error for random intercepts
              , 'eta' #individual random intercept in detection
              , 'eps' #random intercept for site
              , 'psi' #prob augmented individual belongs in N 
              , 'p' #detection probability
              , 'N.site' #estimates of abundance by site
              , 'N.tot' #total abundance across sites and years
)

#initial values for whether an augmented individual is added in
zst <- c( rep(1, N), rep(0, M-N) )
# We select the abundance predictors we want
XIK <- ik_sc[,c("herbaceous")]
#how many ecological predictors
B <- 1#dim(XIK)[2]
#how many detection predictors
A <- 3

#define initial parameter values
inits <- function(){ list( beta = rnorm( B ),
                           alpha = rnorm( A ),
                           int.lam = runif(1), int.p = runif(1),
                           z = zst ) }

#define data that will go in the model
str( win.data <- list( y = as.matrix(y[ ,c("j_1", "j_2", "j_3") ] ),
                       #number of sitesXyear, surveys, det predictors, and abund preds
                       I = I, J = J, A = A,  N = N, M = M,B = B,
                       #siteXyear id for each individual
                       site = y$idno,
                       #site level habitat predictors
                       XIK = XIK,
                       #observation predictors:
                       wind = ij_wide[,widx],
                       temp = ij_wide[,tidx],
                       effort = ij_wide[,eidx]
) )

#call JAGS and summarize posteriors:
ma1 <- autojags( win.data, inits = inits, params, modelname, #
                  n.chains = nc, n.thin = nt, n.burnin = nb,
                  iter.increment = ni, max.iter = 1000000, 
                  Rhat.limit = 1.1,
                  save.all.iter = FALSE, parallel = TRUE ) 

# ma1 <- update( ma1, parameters.to.save= params,
#                    n.iter = 500000, n.thin = nt)

plot(ma1)
summary(ma1)
#####################################################
##### save relevant stuff ##################################
save.image( "MarkedAugResultsJAGs.RData")

################### end of script #######################################