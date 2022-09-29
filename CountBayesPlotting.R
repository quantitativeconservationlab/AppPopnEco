#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
# We import results for Bayesian analysis of single season counts   ###
# which were modelled using a N-mixture model in a Bayesian framework #
#                                                                     #
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
load( "CountBayesResults.RData" )

################################################################################
#################### viewing model results ####################################
##################################################################################
#define model results to plot
mr <- m1

############## trace plots ############
#plot( mr ) #plots traces and posterior densities for all parameters
par( mfrow = c( 2, 2 ), ask = F, mar = c(3,4,2,2) )
#detection parameters
#intercept 
traceplot( mr, parameters = c( 'int.det') )
#random intercepts for observer
traceplot( mr, parameters = c( 'eps.det') )
#fixed effect
traceplot( mr, parameters = c( 'alpha') )
#for abundance
traceplot( mr, parameters = c( 'int.lam') )
traceplot( mr, parameters = c( 'beta') )

############## whisker plots #############
par( mfrow = c( 1,1 ), ask = F , mar = c(3,4,2,2) )
#for detection
#fixed effects for detection
whiskerplot( mr, parameters = c( 'int.det', 'alpha' ) , zeroline = TRUE )
#random intercept for detection
whiskerplot( mr, parameters = c( "eps.det" ) )
#for abundance
par( mfrow = c( 1,1 ), ask = F , mar = c(3,4,2,2) )
#fixed effects for abundance
whiskerplot( mr, parameters = c( 'int.lam', 'beta' ) , zeroline = TRUE )

#derived parameters
whiskerplot( mr, parameters = c( "p" ) )
whiskerplot( mr, parameters = c( "N" ) )


#############################################################
### model evaluation ##########################################

#combine fit statistics into dataframe for plotting
fitdf <- data.frame( fit = mr$sims.list$fit, 
                fit.new = mr$sims.list$fit.new )
head( fitdf)
#plot
ggplot( fitdf, aes( x = fit, y= fit.new ) ) +
  geom_point( size = 2 ) +
  theme_classic( base_size = 17 ) +
  xlab( "Discrepancy actual data" ) +
  ylab( "Discrepancy predicted data" ) +
  geom_abline(intercept = 0, slope = 1)

#calculate fit ratio
mean( mr$mean$fit ) / mean( mr$mean$fit.new )

######### end model evaluation ######################################
##################################################################
######### violin plots ###########################################
# We may want to plot our model coefficient distributions #
# one approach is to use violin plots
# here we demonstrate for the abundance submodel
#start with extracting relevant fixed effects for abundance submodel
beta.mat <- mr$sims.list$beta
#label columns based on order you included predictors in the model
colnames( beta.mat ) <- c( "Cheatgrass (%)", "Sagebrush (%)" )

# convert matrix to long format
beta.df <- tidyr::gather( data = as.data.frame( beta.mat), key = Predictor,
                          Value )
#plot 
ggplot( beta.df, aes( y = Value, x = Predictor, fill = Predictor ) ) + 
  #set plotting theme 
  theme_classic( base_size = 17 ) + 
  #add label for y axis
  ylab( 'Standardized effect size' ) +
  #plot distribution
  geom_violin( trim = FALSE, #fill = '#A4A4A4', 
               size = 0.5 ) +
  #add boxplot highlighting mean and 95 CIs 
  geom_boxplot( width = 0.1 ) +
  #add zero line
  geom_hline( yintercept = 0, size = 1.2, color = 'black' )

#repeat process for detection submodel here:
#
# 
#
####################################################################
######### partial prediction plots #############################
# Estimate partial prediction plots (marginal effect plots) for predictors 
# with 95% CIs not overlapping zero:
# For abundance submodel first:####
# Start by creating our datasets to predict over
# how many values do we use:
n <- 100
#define a vector of ones for intercept
int <- rep( 1, n )
# Use the observed values to define range of predictor:
cheatgrass <- seq( min( closeddf[,"cheatgrass"]),max( closeddf[,"cheatgrass"]),
                   length.out = n )
#standardize predictors:
cheat.std <- scale( cheatgrass )

#extract relevant fixed coefficient from abundance submodel results
fixedabund <- cbind( mr$sims.list$int.lam, mr$sims.list$beta[,1] )

#estimate predicted abundance 
predabund <- exp( fixedabund %*% t( cbind( int, cheat.std) ) )
#calculate mean abundance
mabund <- apply( predabund, MARGIN = 2, FUN = mean )
#calculate 95% Credible intervals for abundance
CIabund <- apply( predabund, MARGIN = 2, FUN = quantile, 
                  probs = c(0.025, 0.975) )

#create dataframe combining all predicted values for plotting
abunddf <- data.frame( mabund, t(CIabund),
                       cheat.std, cheatgrass )
#view
head( abunddf); dim( abunddf)
#rename predicted abundance columns
colnames(abunddf )[1:3] <- c(  "Mean", "lowCI", "highCI" )

#plot marginalised effects for abundance submodel 
ggplot( abunddf, aes( x = cheatgrass, y = Mean) ) +
  theme_classic( base_size = 17) +
  ylab( "Relative abundance" ) +
  xlab( "Cheatgrass (%)" ) +
  geom_line( size = 1.5) +
  geom_ribbon( alpha = 0.3, aes( ymin = lowCI, ymax = highCI ) )

#How does this plot compare to the one you obtained from the 
# unmarked analysis?
# Answer: 
# 
##### detection marginal effects ######
# our only fixed predictor in detection submodel was time
# what are the min max times:
closeddf %>% select( time.j1, time.j2, time.j3 ) %>% 
  summarise_all(list(min, max))
#use them to define your bounds:
Time <- round(seq( 0, 360, length.out = n ),0)
time.std <- scale( Time )
#extract relevant fixed coefficient for detection submodel results
fixeddet <- cbind( mr$sims.list$int.det, mr$sims.list$alpha )
#estimate predicted detection
preddet <- plogis( fixeddet %*% t( cbind( int, time.std ) ) )
#calculate mean abundance
mdet <- apply( preddet, MARGIN = 2, FUN = mean )
#calculate 95% Credible intervals for abundance
CIdet <- apply( preddet, MARGIN = 2, FUN = quantile, 
                  probs = c(0.025, 0.975) )

#create dataframe combining all predicted values for plotting
detdf <- data.frame( mdet, t(CIdet),
                       time.std, Time )
#view
head( detdf); dim( detdf)
#rename predicted abundance columns
colnames(detdf )[1:3] <- c(  "Mean", "lowCI", "highCI" )

#plot marginalised effects for abundance submodel 
ggplot( detdf, aes( x = Time, y = Mean) ) +
  theme_classic( base_size = 17) +
  ylab( "Probability of detection" ) +
  xlab( "Time (mins pass 6:00am)" ) +
  geom_line( size = 1.5) +
  geom_ribbon( alpha = 0.3, aes( ymin = lowCI, ymax = highCI ) )
#How does this plot compare to the one you obtained from the 
# unmarked analysis?
# Answer: 
# 
######## end of marginal effect plots ###################
####################################################################
########## save relevant figures or data ###########################
### 
# Save your figures using what you have learnt in the past:
### Code here: 
###
####################### end of script #################################