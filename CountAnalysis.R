#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we import our cleaned data for 2009 for our point count       ##
#  observations for Piute ground squirrels at the NCA and run a      ##
## closed population N-mixture analysis. The occupancy model is      #
# hierarchical with : (1) an ecological submodel linking abundance to #
## environmental predictors at each site; (2) an observation submodel ##
## linking our detection probability to relevant predictors.         ##
##                                                                   ##
# Female Piute ground squirrels give birth to an average of 5-10 young#
# Reproduction and survival are likely influenced by colder temperatures #
# in Feb, when they come out of hibernation.                          #                
# Survival is also likely affected by hot temperatures, with individuals#
# unable to forage when temperatures are too hot. So we expect a      #
# negative relationship between survival and max T in Apr-May         #
# Survival is expected to be higher in sites with more sagebrush      #
#                                                                     #
# Detection may be related to observer effects and to time of day as a #
# quadratic, with higher detection expected in the middle of the day, #
# when squirrels are expected to be most active.                      #
#######################################################################

##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

#install relevant packages
install.packages( 'Rtools' )
install.packages( "nmixgof" ) #for evaluating N-mixture models

library( tidyverse )#includes dplyr, tidyr and ggplot2
library( unmarked ) #
library( MuMIn )
library( AICcmodavg)
library( nmixgof )
## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Data/", sep = "" )
# load our cleaned data
closeddf <- read.csv( file = paste( datadir, "closed_counts.csv", sep = ""),
                      header = TRUE )
#view
head( closeddf ); dim( closeddf ) 
#### End of data load -------------
####################################################################
##### Ready data for analysis --------------
# What predictors do we think drive colonization, extinction and # 
# detection of Piute ground squirrels at the NCA? #
# Let's define our unmarked dataframe:
# Start by defining which columns represent the response (counts):
umf <- unmarkedFramePCount( y = as.matrix( closeddf[ ,c("count.j1", "count.j2", "count.j3")]),
                          # Define predictors at the site level:
                          siteCovs = closeddf[ ,c("sagebrush", "cheatgrass")],
                          # Define predictors at the survey level as a list:
                          obsCovs = list( obsv = closeddf[ ,c("observer.j1", "observer.j2", "observer.j3")],
                                          time = closeddf[ ,c('time.j1', 'time.j2', 'time.j3') ] ) ) 

# View summary of unmarked dataframe:
summary( umf )
#now scale ecological predictors:
sc <- apply( siteCovs(umf), MARGIN = 2, FUN = scale )
# We replace the predictors in our unmarked dataframe with the scaled values:
siteCovs( umf ) <- sc
#now for observation-level predictors:
osc <- as.vector(scale( obsCovs(umf)[2] ))
#replace with scaled values:
obsCovs(umf)[2] <- osc
#recheck
summary( umf )
### end data prep -----------
### Analyze data ------------------------------------------
# We are now ready to perform our analysis:
fm.closed <- pcount( ~ 1 + obsv + time
                   ~ 1 + sagebrush + cheatgrass, 
                   #we need to define the maximum possible abundance
                   #during the primary occasion
                    K = 1000,
                   data = umf )
# Note that we start with the observation submodel #
#We then define the ecological submodel as related #
# to sagebrush and cheatgrass. We end by defining the data to be used.

# View model results:
fm.closed

# We can also estimate confidence intervals for coefficients in #
# ecological submodel:
confint( fm.closed, type = "state" )
# Why do we call them coefficients and not predictors?
# Answer:
#
# coefficients for detection submodel:
confint( fm.closed, type = 'det' )
#
# Based on the overlap of the 95% CIs for your predictor coefficients, #
# can you suggest which may be important to each of your responses? #
# Answer:
# 
#############end full model ###########
##########################################################################
# Model fit and evaluation -----------------------------------------------

# Now that we looked at the initial output we can evaluate our model to #
# decide if we are happy to proceed or need to modify our analysis #
# somehow.

# We start with goodness of fit (GoF) outlined by Duarte et al. 2018 #
# Ecological modelling 374:51-59 and available via AICmodavg package #
# The Nmix.gof.test relies on a Pearson chi-square to assess the fit of #
# N-mixture models. The approach uses bootstrapping to estimate the p values #
# The test also estimates a c-hat measure of overdispersion, as the  #
# observed test statistic divided by the mean of the simulated test statistics #

# Let's compute observed chi-square, assess significance, and estimate c-hat
gof.boot <- Nmix.gof.test( fm.closed, nsim = 1000, print.table = TRUE )
#view
gof.boot
# What does the output tell us about our model fit?
# Answer:
#
# We also evaluate how well our full model did against the null model # 
# by estimating pseudo-R^2, based on Nagelkerke, N.J.D. (2004) A Note #
# on a General Definition of the Coefficient of Determination. Biometrika 78,#
# pp. 691-692.#
# We run the null model
fm.null <- pcount( ~ 1 ~ 1,
            K = 1000, data = umf )
#view
fm.null
# Now build the list with the two models of interest:
rms <- fitList( 'full' = fm.closed,
                'null' = fm.null )
# Then use model selection function from unmarked, defining which is the null:
unmarked::modSel(rms, nullmod = "null" )


# Now we used the gof checks outlined in Knape et al. 2018 MEE 9:2102-2114
# We start by estimating overdispersion metrics 
chat( fm.closed, type = 'marginal' )
chat( fm.closed, type = 'site-sum' )
chat( fm.closed, type = 'observation' )

# Plot rq residuals against untransformed numeric predictors. This may help
# detect problems with the functional form assumed in a model
residcov( fm.closed )
# What do the plots tell you?
# Answer:
#
# Plot residuals against fitted values. Site-sum randomized quantile residuals
# are used for site covariates while marginal residuals are used for
# observation covariates. 
residfit( fm.closed, type = 'site-sum' )
# Plot the observation model residuals
residfit( fm.closed, type = 'observation' )
# What did Knape et al. 2018 say these residuals were useful for?
# Answer:
#
# Qq plots of randomized residuals against standard normal quantiles #
# Under a good fit residuals should be close to the identity line. 
residqq( fm.closed, type = 'site-sum' )
residqq( fm.closed, type = 'observation' )

# What do these plots indicate? 
# Answer:
# 
# What is some of the advice recommended by Knape et al. 2018 if we want 
# to use abundance estimates from these N-mixture models?
# Answer:
#
# Now try fitting other functional forms (e.g. ZIP or Negative Binomial)
# Do you get a better fit?
# Answer:
#
#########################################################################
##### Summarizing model output ##############

# We now see the effects that our predictors are having on this trend. #
# by plotting partial prediction plots for our ecological submodels #
# Here I focus only on those with 95% CIs not overlapping zero:
# We start by creating our datasets to predict over
# how many values do we use:
n <- 100
# we use the observed values to define our range:
cheatgrass <- seq( min( closeddf[,"cheatgrass"]),max( closeddf[,"cheatgrass"]),
                   length.out = n )
# what are the min max times:
closeddf %>% select( time.j1, time.j2, time.j3 ) %>% 
  summarise_all(list(min, max)) 
Time <- round(seq( 0, 360, length.out = n ),0)
#standardize them
cheat.std <- scale( cheatgrass )
time.std <- scale( Time )
#combine standardized predictor into a new dataframe to predict partial relationship
# with sagebrush. We replace value of other predictor with its mean
abundData <- data.frame( sagebrush = 0, cheatgrass = cheat.std )

#predict partial relationship between sagebrush and occupancy
pred.cheat <- predict( fm.closed, type = "state", newdata = abundData, 
                          appendData = TRUE )
#view
head( pred.cheat ); dim( pred.cheat )

# now for detection
detData <- data.frame( obsv = list(0,0,0,0), time = time.std )
#predict partial relationship for detection:
pred.time <- predict( fm.closed, type = "det", newdata = detData, 
                           appendData = TRUE )

# create plots for ecological submodel:
cheatp <- cbind( pred.cheat[,c("Predicted", "lower", "upper") ], cheatgrass ) %>%
  # define x and y values
  ggplot(., aes( x = cheatgrass, y = Predicted ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Cheatgrass (%)", y = "Relative abundance" ) +
  # add band of confidence intervals
  geom_smooth( aes(ymin = lower, ymax = upper ), 
               stat = "identity",
               size = 1.5, alpha = 0.5, color = "grey" ) +
  # add mean line on top
  geom_line( size = 2 ) 
#view
cheatp
# How do you interpret this relationship?
# Answer:
#
############################################################################
################## Save your data and workspace ###################

# Save workspace:
save.image( "CountResults.RData" )

#save the plot objects you need for your presentation
#Cheatgrass x abundance plot:
tiff( 'Data/CheatXAbund.tiff',
      height = 10, width = 12, units = 'cm', 
      compression = "lzw", res = 400 )
#call the plot
cheatp
#end connection
dev.off()

########## End of saving section ##################################

############# END OF SCRIPT #####################################