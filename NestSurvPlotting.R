#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                  ###
## This script takes model results from NestSurvAnalysis.R which    ###
## models chick survival using data from unmarked individuals.      ##
##
# We visualize model results and model fit in this script            #
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
###################################################################
#### Load or create data -----------------------------------------
#load relevant workspace
load( "SurvivalResults.RData" )

################################################################################
#################### viewing model results ####################################
##################################################################################
#define model results to plot
mr <- m1
#Total number of iterations ran:
N <- mr$mcmc.info$n.samples

#extract summary results
summary( mr )
plot( mr )
############## trace plots ############
#plot( mr ) #plots traces and posterior densities for all parameters
par( mfrow = c( 2, 2 ), ask = F, mar = c(3,4,2,2) )
#detection parameters
#intercept 
traceplot( mr, parameters = c( 'int.p') )
#random intercepts
traceplot( mr, parameters = c( 'sigma.eps.p.M') )
traceplot( mr, parameters = c( 'sigma.eps.p.J') )
traceplot( mr, parameters = c( 'eps.p.M') )
traceplot( mr, parameters = c( 'eps.p.J') )

#survival parameters
#intercept
traceplot( mr, parameters = c( 'int.phi') )
#fixed effects
traceplot( mr, parameters = c( 'beta') )

############## whisker plots #############
par( mfrow = c( 1,1 ), ask = F , mar = c(3,4,2,2) )
#for detection
#fixed effects for detection
whiskerplot( mr, parameters = c( 'int.p' ) , zeroline = TRUE )
#random intercept for detection
whiskerplot( mr, parameters = c( "eps.p.M" ) )
whiskerplot( mr, parameters = c( "eps.p.J" ) )
#for survival
par( mfrow = c( 1,1 ), ask = F , mar = c(3,4,2,2) )
#fixed effects for abundance
whiskerplot( mr, parameters = c( 'int.phi', 'beta' ) , 
             zeroline = TRUE )

#survival at 6 weeks
whiskerplot( mr, parameters = c( 'surv6wk') )
##### end whiskeplots #####

###################################################################################
################# plot partial-predicted relationships ###########
############
head( chickdf ); head( XN )

#define range for important variables in survival model:
#hatching date
hatch <- seq( min( chickdf$hatchDay), max( chickdf$hatchDay ), 
              by = 1 )
sclhatch <- scale( hatch )
#chick age:
age <- seq( 1, length( hatch ), by = 1 )
sclage <- scale( age )

#distance to eagle nest:
dis <- seq( min( chickdf$MinDistOccEgl), max( chickdf$MinDistOccEgl), 
            length.out =  length( hatch ) )
scldis <- scale( dis )

#delta: inverse weight of eagle nest density:
delta <- seq( min( chickdf$DeltaOccEgl), max( chickdf$DeltaOccEgl), 
              length.out =  length( hatch ) )
scldelta <- scale( delta )

#intercept for covariate matrix 
phiint <- rep( 1, length( hatch ) )
#join intercept with relevant predictor components
phi.hatch <- cbind( mr$sims.list$int.phi, mr$sims.list$beta[,1] ) %*%
              t( cbind( phiint, sclhatch ) )
#calculate its mean and get inverse logit
phi.hatch.m <- plogis( apply( phi.hatch, MARGIN = 2, FUN = mean ) )
#calculate 95% CIs and get inverse logit:
phi.hatch.CI <- plogis( apply( phi.hatch, MARGIN = 2,
                   FUN = quantile, probs = c(0.025, 0.975) ) )
# for age
phi.age <- cbind( mr$sims.list$int.phi, mr$sims.list$beta[,4] ) %*%
            t( cbind( phiint, sclage ) )
#calculate its mean and get inverse logit
phi.age.m <- plogis( apply( phi.age, MARGIN = 2, FUN = mean ) )
#calculate its 95% CIs and get its inverse logit:
phi.age.CI <- plogis( apply( phi.age, MARGIN = 2,
               FUN = quantile, probs = c(0.025, 0.975) ) )

# for distance to eagle nest
phi.dis <- cbind( mr$sims.list$int.phi, mr$sims.list$beta[,2] ) %*%
            t( cbind( phiint, scldis ) )
#calculate its mean and get inverse logit
phi.dis.m <- plogis( apply( phi.dis, MARGIN = 2, FUN = mean ) )
#calculate its 95% CIs and get its inverse logit:
phi.dis.CI <- plogis( apply( phi.dis, MARGIN = 2,
               FUN = quantile, probs = c(0.025, 0.975) ) )
# for delta of eagle nest
phi.delta <- cbind( mr$sims.list$int.phi, mr$sims.list$beta[,3] ) %*%
  t( cbind( phiint, scldelta ) )
#calculate its mean and get inverse logit
phi.delta.m <- plogis( apply( phi.delta, MARGIN = 2, FUN = mean ) )
#calculate its 95% CIs and get its inverse logit:
phi.delta.CI <- plogis( apply( phi.delta, MARGIN = 2,
              FUN = quantile, probs = c(0.025, 0.975) ) )

#combine predicted estimates into a dataframe
phi.rships <- data.frame( sclhatch, hatch, sclage, age, 
                          scldis, dis, scldelta, delta,
                          phi.hatch.m, t( phi.hatch.CI ),
                          phi.age.m, t( phi.age.CI ),
                          phi.dis.m, t( phi.dis.CI ),
                          phi.delta.m, t( phi.delta.CI ) )

#view
head( phi.rships )

#plot 
phip <- ggplot( data = phi.rships ) + theme_classic() +
  theme( legend.position = "none", 
         text = element_text( size = 18 ), 
         axis.line = element_line( size = 1.3 ) ) 
phip

#hatching plot
hatchp <- phip + #+ ylim( 0, 1 ) 
  xlab( "Julian hatching day" ) + #xlim( 0,42 ) +
  geom_line( aes( x = hatch, y = phi.hatch.m ), size = 1.2 ) +
  geom_ribbon( alpha = 0.4, aes( x = hatch, ymin = X2.5., ymax = X97.5. ) ) + #, 
  ylab( "Daily survival probability" ) #+ ylim( 0, 0 ) 

hatchp
#age
agep <- phip + #+ ylim( 0, 1 ) 
  xlab( "Age (days)" ) + xlim( 0,42 ) +
  geom_line( aes( x = age, y = phi.age.m ), size = 1.2 ) +
  geom_ribbon( alpha = 0.4, aes( x = age, ymin = X2.5..1, ymax = X97.5..1 ) ) + #, 
  ylab( "Daily survival probability" ) #+ ylim( 0, 0 ) 
agep
#distance
disp <- phip + #+ ylim( 0, 1 ) 
  xlab( "Distance to eagle nest (km)" ) + #xlim( 0,42 ) +
  geom_line( aes( x = dis, y = phi.dis.m ), size = 1.2 ) +
  geom_ribbon( alpha = 0.4, aes( x = dis, ymin = X2.5..2, ymax = X97.5..2 ) ) + #, 
  ylab( "Daily survival probability" ) #+ ylim( 0, 0 ) 

#delta
deltap <- phip + #+ ylim( 0, 1 ) 
  xlab( "Nearby density of eagle nests" ) + #xlim( 0,42 ) +
  geom_line( aes( x = delta, y = phi.delta.m ), size = 1.2 ) +
  geom_ribbon( alpha = 0.4, aes( x = delta, ymin = X2.5..3, ymax = X97.5..3 ) ) + #, 
  ylab( "Daily survival probability" ) #+ ylim( 0, 0 ) 

###### end of pred rship plots #####
