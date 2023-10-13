#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
# We import results for Bayesian analysis of multi-season counts   ###
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
load( "RobCountBayesResults.RData" )

################################################################################
#################### viewing model results ####################################
##################################################################################

#define model results to plot
mr <- m2

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
traceplot( mr, parameters = c( 'gamma') )

############## whisker plots #############
par( mfrow = c( 1,1 ), ask = F , mar = c(3,4,2,2) )
#for detection
#fixed effects for detection
whiskerplot( mr, parameters = c( 'int.det', 'alpha' ) , zeroline = TRUE )
#random intercept for observer effect
whiskerplot( mr, parameters = c( "eps.det" ) )
#for abundance
par( mfrow = c( 1,1 ), ask = F , mar = c(3,4,2,2) )
#fixed effects for abundance
whiskerplot( mr, parameters = c( 'beta' ) , zeroline = TRUE )

#derived parameters
whiskerplot( mr, parameters = c( "p" ) )
whiskerplot( mr, parameters = c( "N" ) )

##################################################################
######### annual plots  ###########################################
# We look at annual changes in abundance for our sites ###

Ndf <- opendf %>% 
  dplyr::select( o.sites, year, siteid, yearid )
#view
head( Ndf ); dim( Ndf )
# start by extracting mean Abundance for your sites
Nwide <- as.data.frame(mr$mean$N)
colnames( Nwide) <-  sort(unique(opendf$year)) 
head(Nwide)
Nwide$siteid <- sort(unique(opendf$siteid ) )
Nlong <- pivot_longer( Nwide, cols = starts_with("2"),
              names_to = "year",
              values_to = "N" )
head(Nlong)
#now add credible intervals
L <- H <- data.frame( mr$mean$N )
for(i in 1:I){
  for(k in 1:K) {
    L[i,k] <-  quantile(
      mr$sims.list$N[ ,i,k ],
      probs = 0.025, na.rm = TRUE )
    H[i,k] <-  quantile(
      mr$sims.list$N[ ,i,k ],
      probs = 0.975, na.rm = TRUE )
  }}
colnames( L ) <- colnames( H ) <-  sort(unique(opendf$year)) 
head(L)
low <- L %>%
  pivot_longer( cols = starts_with("2"),
                names_to = "year",
                values_to = "Nlow" )

high <- H %>%
  pivot_longer( cols = starts_with("2"),
                names_to = "year",
                values_to = "Nhigh" )
Nlong$Nlow <- low$Nlow
Nlong$Nhigh <- high$Nhigh

bind_cols( Nlong, opendf[ ,c("count.j1", "count.j2", "count.j3")] )
#look at annual changes in abundance
ggplot( Nlong, aes( x = year, y = N, 
  group = as.factor(siteid), color= as.factor(siteid)))  +
  theme_bw( base_size = 15 ) +
  theme( legend.position = "none" ) +
  geom_point( size = 3 ) +
  geom_line( linewidth = 2 )

#compare to max counts
ggplot( Nlong, aes( x = year, y = N, 
                  group = as.factor(siteid), 
                  color= as.factor(siteid)))  +
  theme_bw( base_size = 15 ) +
  theme( legend.position = "none" ) +
  geom_point( size = 3 ) +
   geom_line( linewidth = 2 ) +
   geom_ribbon( alpha = 0.3, 
       aes( ymin = Nlow, ymax = Nhigh,
            fill = as.factor(siteid) ) )

opendf %>% 
  dplyr::select(siteid, year, 
                count.j1, count.j2, count.j3) %>% 
  mutate( year = as.character(year) ) %>% 
    left_join( Nlong, by = c("siteid", "year") )

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
cheatgrass <- seq( min( opendf[,"cheatgrass"]),max( opendf[,"cheatgrass"]),
                   length.out = n )
#standardize predictors:
cheat.std <- scale( cheatgrass )

#extract relevant fixed coefficient from abundance submodel results
fixedabund <- cbind( mr$sims.list$int.lam, mr$sims.list$beta[,1] )

#estimate predicted abundance 
predabund <- exp( (fixedabund %*% t( cbind( int, cheat.std ) ) ) +
                    mean(mr$mean$N ) * mr$mean$gamma )
            
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

#repeat for AprMay temp
# Use the observed values to define range of predictor:
aprmayT <- seq( min( opendf[,"AprMay.maxT"]),max( opendf[,"AprMay.maxT"]),
                   length.out = n )
#standardize predictors:
aprmayT.std <- scale( aprmayT )

#extract relevant fixed coefficient from abundance submodel results
aprmayTabund <- cbind( mr$sims.list$int.lam, mr$sims.list$beta[,4] )

#estimate predicted abundance 
tabund <- exp( aprmayTabund %*% t( cbind( int, aprmayT.std) ) +
                 mean(mr$mean$N ) * mr$mean$gamma )
#calculate mean abundance
tmabund <- apply( tabund, MARGIN = 2, FUN = mean )
#calculate 95% Credible intervals for abundance
tCIabund <- apply( tabund, MARGIN = 2, FUN = quantile, 
                  probs = c(0.025, 0.975) )

#create dataframe combining all predicted values for plotting
tabunddf <- data.frame( tmabund, t(tCIabund),
                       aprmayT.std, aprmayT )
#view
head( tabunddf); dim( tabunddf)
#rename predicted abundance columns
colnames(tabunddf )[1:3] <- c(  "Mean", "lowCI", "highCI" )

#plot marginalised effects for abundance submodel 
ggplot( tabunddf, aes( x = aprmayT, y = Mean) ) +
  theme_classic( base_size = 17) +
  ylab( "Relative abundance" ) +
  xlab( "Maximum Temperature (Apr-May)" ) +
  geom_line( size = 1.5) +
  geom_ribbon( alpha = 0.3, aes( ymin = lowCI, ymax = highCI ) )


##### detection marginal effects ######
# our only fixed predictor in detection submodel was time
# what are the min max times:
closeddf %>% select( time.j1, time.j2, time.j3 ) %>% 
  summarise_all(list(min, max))
#use them to define your bounds:
Time <- round(seq( 0, 360, length.out = n ),0)
time.std <- scale( Time )
time2.std <- scale( Time^2 )
#extract relevant fixed coefficient for detection submodel results
fixeddet <- cbind( mr$sims.list$int.det, mr$sims.list$alpha )
#estimate predicted detection
preddet <- plogis( fixeddet %*% t( cbind( int, time.std, time2.std ) ) )
#calculate mean abundance
mdet <- apply( preddet, MARGIN = 2, FUN = mean )
#calculate 95% Credible intervals for abundance
CIdet <- apply( preddet, MARGIN = 2, FUN = quantile, 
                probs = c(0.025, 0.975) )

#create dataframe combining all predicted values for plotting
detdf <- data.frame( mdet, t(CIdet),
                     time.std, time2.std, Time )
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