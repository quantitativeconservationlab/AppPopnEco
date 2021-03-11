#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we import our cleaned data for the time-series of point count #
#  observations for Piute ground squirrels at the NCA and run        ##
## robust population N-mixture analyses. The models are hierarchical  #
#  with : (1) an ecological submodel linking abundance to             #
## environmental predictors; (2) an observation submodel linking     ##
##  detection probability to relevant predictors.                    ##
##                                                                   ##
# Female Piute ground squirrels give birth to an average of 5-10 young#
# Reproduction and survival are likely influenced by colder temperature #
# in Feb, when they come out of hibernation.                          #                
# Survival is likely affected by hot temperatures, with individuals   #
# unable to forage when temperatures are too hot.                     #
# Survival is expected to be higher in sites with more sagebrush.     #
#                                                                     #
# Detection may be related to observer effects and to time of day as a #
# quadratic, with higher detection expected in the middle of the day, #
# when squirrels are most active.                                     #
#                                                                    ##
#######################################################################
##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

# load packages:
library( tidyverse )#includes dplyr, tidyr and ggplot2
library( unmarked ) #analyze count and occupancy data
library( AICcmodavg) #gof tests (Duarte et al. 2018)
library( nmixgof ) #more gof tests (Knape et al. 2018)
## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Data/", sep = "" )
# load our cleaned data
opendf <- read.csv( file = paste( datadir, "open_counts.csv", sep = ""),
                      header = TRUE )
#view
head( opendf ); dim( opendf ) 
#### End of data load -------------
####################################################################
##### Ready data for analysis --------------
#create your response dataframe:
yy <- opendf %>%  
  select( o.sites, year, count.j1, count.j2, count.j3 ) %>% 
  pivot_wider( names_from = year, #which criteria is getting turned to wide
               values_from = c("count.j1", "count.j2", "count.j3"), #what columns do we include
               #this rearranges the naming of the columns to use year first 
               names_glue = "{year}_{.value}" )
#check did we get the right number of dimensions?
head( yy ); dim( yy )
#now we sort columns by year instead of survey
yy <- yy %>% dplyr::select( str_sort( colnames(yy), numeric = TRUE ) ) %>%
  # we remove site id
  select( -o.sites )

#now create siteXyear dataframes:
#for sagebrush:
sagebrush <- opendf %>% select( o.sites, year, sagebrush ) %>% 
          spread( key = year, value = sagebrush ) %>% 
  select( -o.sites )
#check
head(sagebrush); dim(sagebrush)
#for Feb temperature:
Feb.minT <- opendf %>% select( o.sites, year, Feb.minT ) %>% 
  spread( key = year, value = Feb.minT ) %>% 
  select( -o.sites )
#check
head(Feb.minT); dim(Feb.minT)
#for Apr-May temperature:
AprMay.maxT <- opendf %>% select( o.sites, year, AprMay.maxT ) %>% 
  spread( key = year, value = AprMay.maxT ) %>% 
  select( -o.sites )
#check
head(AprMay.maxT); dim(AprMay.maxT)

#create observer level dataframes
#starting with observer effect:
observer <- opendf %>%  
  select( o.sites, year, observer.j1, observer.j2, observer.j3 ) %>% 
  pivot_wider( names_from = year, #which criteria is getting turned to wide
               values_from = c("observer.j1", "observer.j2", "observer.j3"), #what columns do we include
               #this rearranges the naming of the columns to use year first 
               names_glue = "{year}_{.value}" )
#check did we get the right number of dimensions?
head( observer ); dim( observer )
#now we sort columns by year instead of survey
observer <- observer %>% dplyr::select( str_sort( colnames(observer), numeric = TRUE ) ) %>%
  # we remove site id
  select( -o.sites )
#now for time of day:
Time <- opendf %>%  
  select( o.sites, year, time.j1, time.j2, time.j3 ) %>% 
  pivot_wider( names_from = year, #which criteria is getting turned to wide
               values_from = c("time.j1", "time.j2", "time.j3"), #what columns do we include
               #this rearranges the naming of the columns to use year first 
               names_glue = "{year}_{.value}" )
#check did we get the right number of dimensions?
head( Time ); dim( Time )
#now we sort columns by year instead of survey
Time <- Time %>% dplyr::select( str_sort( colnames(Time), numeric = TRUE ) ) %>%
  # we remove site id
  select( -o.sites )

#lastly calculate number of primary periods
T <- length( unique( opendf$year ) )
#now combine into unmarked dataframe:
umf <- unmarkedFramePCO( y = as.matrix(yy),
                        #combine yearXsite dataframes into a list
                        yearlySiteCovs = list( sagebrush = sagebrush,
                                               Feb.minT =Feb.minT,
                                               AprMay.maxT =  AprMay.maxT),
                        #combine observer-level dataframes into list
                        obsCovs = list( observer = observer, 
                                        time = Time ),
                        #provide number of primary periods
                        numPrimary = T )
# Check
summary( umf )
#now scale ecological predictors:
sc <- apply( yearlySiteCovs(umf), MARGIN = 2, FUN = scale )
# We replace the predictors in our unmarked dataframe with the scaled values:
yearlySiteCovs( umf ) <- sc
#now for observation-level predictors:
osc <- as.vector(scale( obsCovs(umf)[2] ))
#replace with scaled values:
obsCovs(umf)[2] <- osc
#recheck
summary( umf )
### end data prep -----------
########################################################################
### Analyze data ------------------------------------------
# The pcountOpen() function allows a wide variety of model parametizations. #
# With more choices, comes potential confusion as to how to implement different #
# options adequately. Models are also very computationally-intensive and can #
# take several hours (or days) to run, depending on your computer capabilities. # 
# If you can provide initial values, do so. Minimize your K in the first run and#
#  start without estimating standard errors. #
# This will allow you to check that you have specified the model correctly. #
# Then build up from there. 

# Here we start with description of some of the options available. Please #
# check the unmarked manual for further details. #

# Initial abundance is always: N[i,t] ~ Poisson(lambda) #
# Survivors can be estimated from a Binomial process as density-dependent:
# S[i,t+1] ~ Binomial( N[i,t], omega )
# if dynamics = "autoreg":
# Recruits can be related to previous abundance, or density-dependent as:
# R[i,t+1] ~ Poisson( N[i,t] * gamma ), where gamma represents per-capita recruitment
# OR:
# if dynamics = 'constant':
# R[i,t+1] ~ Poisson( gamma ), where gamma becomes constant recruitment rate
# Note that under this second option, extinct sites (N[i,t]=0) cannot experience 
# recruitment (i.e, cannot be rescued).
#
# If dynamics = "notrend": lambda * (1-omega) and there is no temporal trend
# If dynamics = "trend": then population growth is models as exponential growth:
# N[i,t] = N[i,t-1] * gamma, where gamma is finite rate of increase (normally #
# referred to as lambda ).
# Dynamics can also be modeled as density-dependent in two other ways:
# where gamma is maximum instantaneous population growth rate (r) and 
# omega is equilibrium abundance (normally termed K).
# (1) Ricker where N[i,t] = N[i,t-1] * exp( gammma * (1 - N[i,t-1] / omega ) ), OR
# (2) Gompertz where N[i,t] = N[i,t-1] * exp( gamma * ( 1 - log( N[i,t-1] + 1 ) /
#     log( omega + 1) ) )

# Let's start with a simple model
fm.notrend <- pcountOpen( #lambda formula for initial abundance:
  lambdaformula = ~1, 
  #gamma is modeled as lambda * ( 1 - omega)
  gammaformula = ~1, 
  #omega 
  omegaformula = ~1, 
  #detection formula:
  pformula = ~1, 
  #no trend
  dynamics = 'notrend', 
  #Define the maximum possible abundance
  K = 400,
  # don't calculate standard errors, which makes it run faster:
  se = FALSE, #useful for the first run
  # set distribution as Poisson:
  mixture = "P", #NB or ZIP also possible 
  # set to true if you want to separate migration from survival:
  immigration = FALSE,
  # provide data
  data = umf, 
  # set control parameters so that you can see progress:
  control = list( trace = TRUE, REPORT = 1) )

# View model results:
fm.notrend

# Now let's simulate population growth for the following years using a #
# Gompertz model adapted to discrete time steps. #
# See: Cruz et al. 2013 PLOS ONE 8(9):e73544 for example.
fm.gompertz <- pcountOpen( #lambda formula for initial abundance:
  lambdaformula = ~1, 
  #gamma is instantaneous population growth rate (r)
  gammaformula = ~1 + sagebrush + Feb.minT + AprMay.maxT, 
  #omega is carrying capacity (K)
  omegaformula = ~1, 
  #detection formula:
  pformula = ~1 + observer + time + I(time)^2,  
  #density-dependent population growth: N[i,t] = N[i,t-1] * 
  #exp( gamma * (1 - log(N[i,t-1] + 1) ) / log(omega + 1) )
  dynamics = 'gompertz', 
  #Define the maximum possible abundance
  K = 500,
  # doesn't calculate standard errors, which makes it run faster:
  se = FALSE, #useful for the first run
  # set distribution as Poisson:
  mixture = "P", #NB or ZIP also possible 
  # set to true if you want to separate migration from survival:
  immigration = FALSE,
  # provide data
  data = umf, 
  # set control parameters so that you can see progress:
  control = list( trace = TRUE, REPORT = 1) )

# View model results:
fm.gompertz

#backtransform parameter estimates
lam <- exp(coef(fm.gompertz, type = "lambda"))
om <- plogis(coef(fm.gompertz, type = "omega"))
gam <- exp( coef( fm.gompertz, type = "gamma"))
p <- plogis( coef( fm.gompertz, type = "det" ) )
#view
lam; om; gam; p
# We can also estimate confidence intervals for coefficients in #
# ecological submodel:
confint( fm.gompertz, type = "lambda" )
confint( fm.gompertz, type = "omega" )
confint( fm.gompertz, type = "gamma" )
# coefficients for detection submodel:
confint( fm.gompertz, type = 'det' )
#
# Based on the overlap of the 95% CIs for your predictor coefficients, #
# can you suggest which may be important to each of your responses? #
# Answer:
#

# What would the Dail, Madsen (2011) model look like?
# Answer:
#
#############end full model ###########
##########################################################################
# Model fit and evaluation -----------------------------------------------
# We start with goodness of fit (GoF) outlined by Duarte et al. 2018 #
# Ecological modelling 374:51-59 and available via AICmodavg package #
# The Nmix.gof.test relies on a Pearson chi-square to assess the fit of #
# N-mixture models. The approach uses bootstrapping to estimate the p values #
# The test also estimates a c-hat measure of overdispersion, as the  #
# observed test statistic divided by the mean of the simulated test statistics #

# Let's compute observed chi-square, assess significance, and estimate c-hat
gof.boot <- Nmix.gof.test( fm.gompertz, nsim = 500, print.table = TRUE )
#view
gof.boot
# Remember that higher chi-squared values represent worse fit

# What about a comparison of our fitted vs observed values?
plot(  yy[,1], fitted( opendf )[,1] )
for( j in 2:dim(y)[2]){
  points( yy[,j], fitted( opendf )[,j] )
}
# Now use gof checks outlined in Knape et al. 2018 MEE 9:2102-2114
# We start by estimating overdispersion metrics 
chat( fm.gompertz, type = 'marginal' )
chat( fm.gompertz, type = 'site-sum' )
chat( fm.gompertz, type = 'observation' )

# Plot rq residuals against untransformed numeric predictors. This may help
# detect problems with the functional form assumed in a model
residcov( fm.gompertz )
# What do the plots tell you?
# Answer:
#
# Plot residuals against fitted values. Site-sum randomized quantile residuals
# are used for site covariates while marginal residuals are used for
# observation covariates. 
residfit( fm.gompertz, type = 'site-sum' )
# Plot the observation model residuals
residfit( fm.gompertz, type = 'observation' )
# What did Knape et al. 2018 say these residuals were useful for?
# Answer:
#
# Niw plot Qq plots of randomized residuals against standard normal quantiles. #
# Under a good fit, residuals should be close to the identity line. 
residqq( fm.gompertz, type = 'site-sum' )
residqq( fm.gompertz, type = 'observation' )

#########################################################################
##### Summarizing model output ##############
# Now we plot the results we are interested in.
# We can see how mean occupancy changed through time by extracting values from
# the projected.mean data table:
fm.gompertz@projected.mean
# select occupied row and combine it with year
data.frame( year = sort( unique( opendf$year) ),
            N = as.vector( fm.gompertz@projected.mean["abundance",] ) ) %>%
  ggplot(., aes( x = year, y = N ) ) +
  theme_bw(base_size = 15 ) + 
  geom_line( size = 2 )

# What is happening to Piute ground squirrels at the NCA?
# Answer:
# 
# We now see the effects that our predictors are having on this trend. #
# by plotting partial prediction plots for our ecological submodels #
# Here I focus only on those with 95% CIs not overlapping zero:
# We start by creating our datasets to predict over
# how many values do we use:
n <- 100
# we use the observed values to define our range:
sage <- seq( min( opendf[,"sagebrush"]),max( opendf[,"sagebrush"]),
                  length.out = n )
#standardize them
sage.std <- scale( sage )
#combine standardized predictor into a new dataframe to predict partial relationship
sageData <- data.frame( sagebrush = sage.std, AprMay.maxT = 0, Feb.minT = 0 )
#predict partial relationship
pred.sage <- predict( fm.gompertz, type = "gamma", newdata = sageData, 
                      appendData = TRUE )
#view
head( pred.sage ); dim( pred.sage )

#Feb.minT
minT <- seq( min( opendf[,"Feb.minT"]),max( opendf[,"Feb.minT"]),
             length.out = n )
#standardize them
minT.std <- scale( minT )
#combine standardized predictor into a new dataframe to predict partial relationship
minData <- data.frame( sagebrush = 0, Feb.minT = minT.std, AprMay.maxT = 0)
#predict partial relationship
pred.minT <- predict( fm.gompertz, type = "gamma", newdata = minData, 
                           appendData = TRUE )

#AprMay.maxT
maxT <- seq( min( opendf[,"AprMay.maxT"]),max( opendf[,"AprMay.maxT"]),
             length.out = n )
#standardize them
maxT.std <- scale( maxT )
#combine standardized predictor into a new dataframe to predict partial relationship
maxData <- data.frame( sagebrush = 0, Feb.minT = 0, AprMay.maxT = maxT.std )
#predict partial relationship
pred.maxT <- predict( fm.gompertz, type = "gamma", newdata = maxData, 
                      appendData = TRUE )

# create plots for ecological submodels
sagep <- cbind( pred.sage[,c("Predicted", "lower", "upper") ], sagebrush ) %>%
  # define x and y values
  ggplot(., aes( x = sagebrush, y = Predicted ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Sagebrush (%)", y = "Instantaneous growth rate" ) +
  # add band of confidence intervals
  geom_smooth( aes(ymin = lower, ymax = upper ), 
               stat = "identity",
               size = 1.5, alpha = 0.5, color = "grey" ) +
  # add mean line on top
  geom_line( size = 2 ) 
#view
sagep
# How do you interpret this relationship?
# Is there a potential threshold beyond which colonization becomes unlikely?
# Answer:
#
minTp <- cbind( pred.minT[,c("Predicted", "lower", "upper") ], minT ) %>%
  # define x and y values
  ggplot(., aes( x = minT, y = Predicted ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Minimum temperature (Feb)", y = "Instantaneous growth rate" ) +
  # add band of confidence intervals
  geom_smooth( aes(ymin = lower, ymax = upper ), 
               stat = "identity",
               size = 1.5, alpha = 0.5, color = "grey" ) +
  # add mean line on top
  geom_line( size = 2 ) 
#view
minTp
maxTp <- cbind( pred.maxT[,c("Predicted", "lower", "upper") ], maxT ) %>%
  # define x and y values
  ggplot(., aes( x = maxT, y = Predicted ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Maximum temperature (Apr-May)", y = "Instantaneous growth rate" ) +
  # add band of confidence intervals
  geom_smooth( aes(ymin = lower, ymax = upper ), 
               stat = "identity",
               size = 1.5, alpha = 0.5, color = "grey" ) +
  # add mean line on top
  geom_line( size = 2 ) 
#view
maxTp
# Add plots for detection
# Answer:
#
############################################################################
################## Save your data and workspace ###################

# Save workspace:
save.image( "RobCountResults.RData" )

#save the plot objects you need for your presentation
#start by calling the file where you will save it
tiff( 'Data/SageXGam.tiff',
      height = 10, width = 12, units = 'cm', 
      compression = "lzw", res = 400 )
#call the plot
sagep
#end connection
dev.off()

########## End of saving section ##################################

############# END OF SCRIPT #####################################