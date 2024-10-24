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
# We will focus on linking predictors to gains instead of survival   #
# Female Piute ground squirrels give birth to an average of 5-10 young#
# Reproduction is likely more successful when the season last longer #
# so warmer Feb temperatures will suggest quicker end of hibernation. #                
# If Apr-May temperatures are too hot it may decrease food for young. #
# Areas with more sagebrush may have more food sources for young.      #
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
# Be on the lookout for strange results such as likelihoods that don't decrease #
# when you add additional parameters, or NaNs when se=TRUE.

# Here we start with description of some of the options available. Please #
# check the unmarked manual for further details. #

# Initial abundance is always: N[1,t] ~ Poisson(lambda) #
# Survivors can be estimated from a Binomial process as density-dependent:
# S[i,t+1] ~ Binomial( N[i,t], omega ), where omega represents apparent survival
# if dynamics = "autoreg":
# Recruits can be related to previous abundance, or density-dependent as:
# R[i,t+1] ~ Poisson( N[i,t] * gamma ), where gamma represents per-capita recruitment
# Note that under this second option, extinct sites (N[i,t]=0) cannot experience 
# recruitment (i.e, cannot be rescued).
# OR:
# if dynamics = 'constant':
# R[i,t+1] ~ Poisson( gamma ), where gamma becomes constant recruitment rate
#
# If dynamics = "notrend": lambda * (1-omega) so there is no gamma or temporal trend 
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
fm.0 <- pcountOpen( #lambda formula for initial abundance:
  lambdaformula = ~1, 
  #gamma is modeled as lambda * ( 1 - omega)
  gammaformula = ~1, 
  #omega 
  omegaformula = ~1, 
  #detection formula:
  pformula = ~1, 
  #no trend
  dynamics = 'autoreg', 
  #Define the maximum possible abundance
  K = 80,
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
fm.0

# We now run a dynamic count model from Dail-Madsen (2011) model. #
fm.dyn0 <- pcountOpen( #lambda formula for abundance in year 1 only:
  lambdaformula = ~1, 
  #gamma is recruitment rate
  gammaformula = ~1, 
  #omega is survival probability 
  omegaformula = ~1,
  #detection formula:
  pformula = ~1,  #
  #trend
  dynamics = 'constant', 
  #upper bound of discrete integration. Should be higher than the maximum 
  #observed count and high enough that it does not affect #
  #parameter estimates 
  K = 80,
  #for the first run turn off calculation of standar errors
  se = FALSE, 
  #more complicated models allow immigration to be split 
  # from births # we don't enable that here
  immigration = FALSE, 
  #here you provide your unmarked dataframe
  data = umf, 
  control = list( trace = TRUE, REPORT = 1) )
#view
fm.dyn0

# We rerun the model above by adding predictors. 
# The model is hard to fit so we start with an intercept only version #
# and we also start without calculating standard errors.#
fm.dyn <- pcountOpen( #lambda formula for initial abundance:
  lambdaformula = ~1, 
  gammaformula = ~1 + sagebrush + Feb.minT + AprMay.maxT, 
  omegaformula = ~1,
  pformula = ~1 + time,
  dynamics = 'constant', 
  K = 80, 
  se = FALSE, 
  immigration = FALSE,
  data = umf, 
  control = list( trace = TRUE, REPORT = 1) )

# View model results:
fm.dyn

# Next we rerun the model with se=TRUE and use the values of the previous #
# run, as initial values for the new model
#define coefficients from the previous model as initial values for the next run
inits <- coef( fm.dyn )
#run model
fm.dynSE <- pcountOpen( lambdaformula = ~1, 
        gammaformula = ~1 + sagebrush + Feb.minT + AprMay.maxT,
        omegaformula = ~1,  pformula = ~1 + time, 
        K = 80, se = TRUE, data = umf,
        control = list( trace = TRUE, REPORT = 1), starts = inits )

#define which model you want to view
fm <- fm.dynSE
#view results
summary(fm)
#backtransform parameter estimates
exp(coef(fm, type = "lambda"))
#depending on your dynamics, the transformation is a expit or exp:
plogis(coef(fm, type = "omega"))
#exp(coef(fm, type = "omega"))
#depending on your dynamics, gamma may not be present (e.g. in no trend model)
exp( coef( fm, type = "gamma"))
plogis( coef( fm, type = "det" ) )

# We can also estimate confidence intervals for coefficients in #
# ecological submodel:
# 1st season mean abundance:
confint( fm, type = "lambda" )
# omega is apparent survival in the no trend model, K in Gompertz
confint( fm, type = "omega" )
# There is no gamma in the notrend model,  
confint( fm, type = "gamma" )
# coefficients for detection submodel:
confint( fm, type = 'det' )
#
# Based on the overlap of the 95% CIs for your predictor coefficients, #
# can you suggest which may be important to each of your responses? #
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
gof.boot <- Nmix.gof.test( fm, nsim = 100, print.table = TRUE )
#view
gof.boot
# Remember that higher chi-squared values represent worse fit

#########################################################################
##### Summarizing model output ##############
# We now check the effects that our predictors are having on abundance #
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
pred.sage <- predict( fm, type = "gamma", newdata = sageData, 
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
pred.minT <- predict( fm, type = "gamma", newdata = minData, 
                           appendData = TRUE )

#AprMay.maxT
maxT <- seq( min( opendf[,"AprMay.maxT"]),max( opendf[,"AprMay.maxT"]),
             length.out = n )
#standardize them
maxT.std <- scale( maxT )
#combine standardized predictor into a new dataframe to predict partial relationship
maxData <- data.frame( sagebrush = 0, Feb.minT = 0, AprMay.maxT = maxT.std )
#predict partial relationship
pred.maxT <- predict( fm, type = "gamma", newdata = maxData, 
                      appendData = TRUE )

# create plots for ecological submodels
sagep <- cbind( pred.sage[,c("Predicted", "lower", "upper") ], sage ) %>%
  # define x and y values
  ggplot(., aes( x = sage, y = Predicted ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Sagebrush (%)", y = "Recruitment" ) +
  # add band of confidence intervals
  geom_smooth( aes(ymin = lower, ymax = upper ), 
               stat = "identity",
               size = 1.5, alpha = 0.5, color = "grey" ) +
  # add mean line on top
  geom_line( size = 2 ) 
#view
sagep
# How do you interpret this relationship?
# Answer:
#
minTp <- cbind( pred.minT[,c("Predicted", "lower", "upper") ], minT ) %>%
  # define x and y values
  ggplot(., aes( x = minT, y = Predicted ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Minimum temperature (Feb)", y = "Recruitment" ) +
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
  labs( x = "Maximum temperature (Apr-May)", y = "Recruitment" ) +
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
# tiff( 'Data/SageXGam.tiff',
#       height = 10, width = 12, units = 'cm', 
#       compression = "lzw", res = 400 )
# #call the plot
# sagep
# #end connection
# dev.off()

########## End of saving section ##################################

############# END OF SCRIPT #####################################