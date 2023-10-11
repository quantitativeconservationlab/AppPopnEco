#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                  ###
### The scripts readies data for analysis of chick survival in a    ###
## Bayesian framework. The analysis was part of Cruz et. al 2023    ## 
## Avian Conservation & Ecology manuscript evaluating the potential ##
## impacts of Bald Eagles on survival of Common Loon chicks once they ##
## enter the water the day of hatching.Chicks were monitored over three ##
### years from when entered the water until they were expected to fletch #
# at 6 weeks of age. When two chicks from a single nest hatched and #
## entered the water successfully, they could not be identified     #
## separately and were assumed to be the same age for a given brood  #

## Data file: Each row displays the status of a chick in a territory, surv = #
#(1 if observed, 0 if not), as well as the estimated hatch date,     #
# age of chick (as number of days), yrid = year id,  # 
# terid = territory ID, id = chick ID, first and last date of monitoring, #
# age = age when it was last detected, MonDay = julian day for first day of #
# monitoring, hatchDay = julian date when it hatched, cummon = survey #
# length (in days) that chick was monitored, number of surveys = visits,#
# j_id = survey id, distance (km) to nearest occupied Bald Eagle nest #
# MinDistOccEgl, and density of surrounding Bald Eagles = DeltaOccEgl. #
#Further detail on methodology, code use for survival analysis and  #
# results are provided in the manuscript.                           #
#######################################################################

##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

# load packages:
library( tidyverse ) 
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( jagsUI )
## end of package load ###############
########### standardising functions  #################
#genetic function to standardise covariates using Gelman's suggestion:
standardise <- function( xmat, stdevs = 2, marg = c( 1, 2, 3 ) ) { 
  mean.xmat = mean( as.vector( xmat ), na.rm = TRUE )
  sd.xmat = sd( as.vector( xmat ), na.rm = TRUE ) 
  std.xmat = apply( xmat, marg, function( x ){
    ( x - mean.xmat ) / (stdevs * sd.xmat ) } )
  return( std.xmat )
}
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are. 
# I have it in a Data folder in my Rstudio project:
datadir <- paste( getwd(), "/Data/Survival/", sep = "" )

# load observed occurrences:
raw_df <- read.csv( file = paste( datadir, "indvdf.csv", sep = ""),
                    header = TRUE )

                                                                                                             ``########## end of data load ################################
#######################################################################
######## prepare data for analysis #############

#create clean dataframe to store data from the imported one
chickdf <- raw_df 
#view
head( chickdf ); dim( chickdf )

##define general parameters 
#number of territories monitored:
M <- max( chickdf$terid )
#year range for longterm data:
yrrange <- sort( unique( chickdf$yrid ) ) 
#total number of years sampled:
K <- max( chickdf$yrid )
# #max number of days a territory was monitored
# J <- max( chickdf$CumMon )
#OR if you want to truncate to a particular period then 
# specify that here
# We want survival for the first 6 weeks
J <- 42
#number of individual chicks
N <- max( chickdf$id )
### end general parameters 

#prepare matrix of observations:
#create an empty 2D dataframe of rows  = N = number of chicks, and columns = #
# J max number of surveys for a territory # 
y_obs <- array( NA, dim = c(  N , J ) )
#define empty predictor dataframes 
#loop through each chick and extract observed survival values for the 
#days when monitoring occurred. 
#first truncate max monitoring to length J
chickdf$CumMon[ which(chickdf$CumMon > J) ] <- J
for( i in 1:dim( chickdf )[1] ){
  y_obs[ as.numeric( chickdf[ i,'id' ] ),
         as.numeric( chickdf[ i,'CumMon' ] ) ] <- chickdf[ i, "surv" ]
}

#check that it worked:
y_obs[1:5,]

#backfill known survival
z_obs <- array( NA, dim = c(  N , J ) )
#loop to fill out known survival up to last observed:
for( n in 1:N ){
  end <- dplyr::last( which( y_obs[ n,  ] == 1 ) )
  if( !is.na( end ) ) { z_obs[ n, 1:end ] <- 1 }
}
#make first survey 1
z_obs[ , 1 ] <- 1

#extract site-level (territory) covariates
head( chickdf )

# for analysis we will focus on individual hatching date, 
# distance to Bald Eagle nest and density of surrounding Eagles
sitedf <- chickdf %>% 
          #group by individual chick
          group_by( id ) %>% 
          #only keep a single record at this level
          slice( n() ) %>% 
          #select predictors of interest
          dplyr::select( hatchDay, MinDistOccEgl, DeltaOccEgl  ) 

#check output
tail(sitedf); dim(sitedf )

#create matrix of scaled predictors
XN <- apply( sitedf[,c("hatchDay", "MinDistOccEgl", "DeltaOccEgl" )], 
             MARGIN = 2, scale )

#view
head(XN)

### create wide matrices of siteXsurvey predictors 
#age for chick on the day 
agemat <- array( NA, dim = c(  N, J ) )

#create an object that contains details for last monitoring day for each chick
lastdf <- chickdf %>% 
        group_by( id ) %>% 
      slice( n() ) %>% 
       select( terid, age, CumMon )
#check
head( lastdf ); dim( lastdf )
#loop over each individual chick
for( n in 1:dim( lastdf )[1] ){
  c <- pull(lastdf[ n, "CumMon" ])
  print( "number of days monitored")
  print( c )
  #extract age and monitoring day when surveys occurred
  agemat[ n, c ] <- as.numeric( lastdf[ n, "age" ] )

  #fill up rest of the matrix up to max day of surveys
  if( c < J ){
    for( j in (c + 1 ):J){
      agemat[ n, j ] <- agemat[ n, j-1 ] + 1
    }}
  #for 
  for( j in c:1 ){
    agemat[ n, j-1 ] <- agemat[ n, j ] - 1
    #monitor progress
    print( agemat[n, j-1] )
  }#,
}

#check 
agemat[ 1:10, ]

#standardize matrices
AgeMat <- standardise( agemat[,], stdevs = 1, marg = c( 1, 2 ) )
#check
AgeMat[1:10,]

########################################################################
############################################################################
####JAGS code for chick survival with unmarked individuals  ####
#### phi ~ int + hatchday + egldist + egldens  + age   ###
#### p ~ int + eps[m] + eps[j] 
#############################################################################
############## Specify model in bugs language:  #####################
sink( "sm1.txt" )
cat( "
     model{
      #Priors linked to survival 
      #define survival intercept as mean probs:
      int.phi <- log( mean.phi / ( 1 - mean.phi ) )
      #mean survival prob for at least one post-fledgling young
      mean.phi ~ dbeta( 4, 4 ) 
     
      #priors for fixed predictors in survival model:
      for( q in 1:Q ){ #loop over number of predictors
        beta[ q ] ~ dnorm( 0, 0.1 ) 
      } #Q

      # priors linked to detection 
      #detection intercept as mean detection probability:
      int.p <- logit( mean.p )
      #mean prob of detecting a post-fledgling young in the water
      mean.p ~ dbeta( 4, 4 ) 
     
      #random territory intercepts
      for ( m in 1:M ){  #loop over species
        eps.p.M[ m ] ~ dnorm( 0, prec.eps.p.M ) T(-7, 7) 
      } #M
     
      #associated precision of random territory intercepts:     
      prec.eps.p.M <- 1/ ( sigma.eps.p.M * sigma.eps.p.M )
      sigma.eps.p.M ~ dt( 0, 2.5, 7 ) T( 0, )
     
      #random day intercepts
      for ( j in 1:(J-1) ){  #loop over monitoring days
        eps.p.J[ j ] ~ dnorm( 0, prec.eps.p.J ) T(-7, 7) 
      } #J
     
      #associated precision of random intercept for day of monitoring:     
      prec.eps.p.J <- 1/ ( sigma.eps.p.J * sigma.eps.p.J )
      sigma.eps.p.J ~ dt( 0, 2.5, 7 ) T( 0, )
     
      #likelihood
      #ecological survival model
      for( n in 1:N ) { #loop over years
        for( j in 2:J ){ #loop over monitoring days
          #latent, true survival
          z[ n, j ] ~ dbern( phi[ n, j-1 ] * z[ n, j-1 ] ) 
     
          logit( phi[ n, j-1 ] ) <- int.phi + 
              #fixed predictors: 
              #hatching date, distance to eagle nest, inv dist weight to eagle nests
              inprod( beta[ 1:(Q-1) ], XN[ n, 1:(Q-1) ] ) + 
               #age of young
              beta[ Q ] * AgeMat[ n, j-1 ] #
        } #J
      } #N
     
      #observation model for dynamic nest occupancy:
      #loop over seasons 
      for ( n in 1:N ) { 
        #loop over monitoring days
        for( j in 2:J ){ 
     
          logit( p[ n, j-1 ] ) <- int.p + 
                        #random intercept for territory
                        eps.p.M[ siteid[ n ] ] + 
                        #random intercept for monitoring day
                         eps.p.J[ j-1 ]
                         
          #observed survival 
          y_obs[ n, j ] ~ dbern( z[ n, j ] * p[ n, j-1 ] ) 
     
        } #close J loop
      } #close N loop
      
     #Estimate derived estimates
     for( n in 1:N ){
     
      #probability of survival to 6 weeks of age
      surv6wk[ n ] <- prod( phi[ n, 1:( J - 1 ) ] )
      
    }#N
     } #model close
     
     ", fill = TRUE )

sink()
################ end of model specification  #####################################
modelname <- "sm1.txt"
#parameters monitored #only keep those relevant for model comparisons (with different variances)
params <- c( 'surv6wk' #6 week post-fledgling survival
             , 'int.phi' #intercept for survival model 
             , 'beta' #fixed coefficients for survival
             , 'int.p' #intercept for detection model 
             , 'eps.p.M' #random territory intercepts 
             , 'eps.p.J' #random monitoring day intercepts 
             , 'sigma.eps.p.M' #std devs for random intercepts
             , 'sigma.eps.p.J' #std devs for random intercepts
             , 'phi' #estimated survival probability [n,j]
             , 'z' #estimated survival[n,j]
             , 'p' #estimated detection [m,j]
)

#define how many predictors we want to include in the survival model
Q <- 4
#define initial values
inits <- function(){ list( beta = rnorm( Q ) )}

#create list of input data 
str( win.data <- list( y_obs = y_obs, #observed survival
                       z = z_obs, #known daily survival 
                       N = N, M = M, J = J, Q = Q,
                       XN = XN, #individual level predictors
                       AgeMat = AgeMat, 
                       siteid = lastdf$terid #territory id for random intercepts
) )     

#call JAGS and summarize posteriors:
m1 <- autojags( win.data, inits = inits, params, modelname, #
                 n.chains = 3, n.thin = 20,  n.burnin = 50000,
                 iter.increment = 20000, max.iter = 500000, 
                 Rhat.limit = 1.02,
                 save.all.iter = FALSE, parallel = TRUE ) 

m1 <- update( m1, parameters.to.save= params,
                  n.iter = 1000000, n.thin = 20 )
            
###### end sm1 ########

###################  Save relevant objects #################

### save workspace so you can load it into analysis :
save.image( "SurvivalResults.RData" )
##########

####################### end of script ########################################