#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we import our cleaned data for a single year of point count  ##
#  observations for Piute ground squirrels at the NCA and run a      ##
## closed population N-mixture analysis. The model is hierarchical    #
#  with : (1) an ecological submodel linking abundance to             #
## environmental predictors at each site; (2) an observation submodel #
## linking detection probability to relevant predictors.             ##
##                                                                   ##
# Abundance is expected to be higher in sites with more sagebrush     #
# and lower in those with more cheatgrass.                            #                                        #
# Detection may be related to observer effects and to time of day     #
#######################################################################

##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

#install relevant packages
install.packages( 'Rtools' )
install.packages( "nmixgof" ) #for evaluating N-mixture models
#load relevant packages
library( tidyverse )#includes dplyr, tidyr and ggplot2
library( unmarked ) #runs N-mixture models
library( MuMIn ) #calculates pseudo-R^2
library( AICcmodavg) #gof tests (Duarte et al. 2018)
library( nmixgof ) #more gof tests (Knape et al. 2018)
## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Data/", sep = "" )
# load cleaned data
closeddf <- read.csv( file = paste( datadir, "closed_counts.csv", sep = ""),
                      header = TRUE )
#view
head( closeddf ); dim( closeddf ) 
#### End of data load -------------
####################################################################
##### Ready data for analysis --------------
# Let's define our unmarked dataframe:
# Start by defining which columns represent the response (counts):
umf <- unmarkedFramePCount( 
  y = as.matrix( closeddf[ ,c("count.j1", "count.j2", "count.j3")]),
  # Define predictors at the site level:
  siteCovs = closeddf[ ,c("sagebrush", "cheatgrass")],
  # Define predictors at the survey level as a list:
  obsCovs = list( obsv = closeddf[ ,c("observer.j1", "observer.j2", "observer.j3")],
  time = closeddf[ ,c('time.j1', 'time.j2', 'time.j3') ],
  time2 = cbind( (closeddf$time.j1)^2, (closeddf$time.j2)^2, 
                                     (closeddf$time.j3)^2 ) ) ) 

# View summary of unmarked dataframe:
summary( umf )
#now scale ecological predictors:
sc <- apply( siteCovs(umf), MARGIN = 2, FUN = scale )
# We replace the predictors in our unmarked dataframe with the scaled values:
siteCovs( umf ) <- sc
#now for time:
timesc <- as.vector(scale( obsCovs(umf)[2] ))
#replace with scaled values:
obsCovs(umf)[2] <- timesc
time2sc <- as.vector(scale( obsCovs(umf)[3] ))
#replace with scaled values:
obsCovs(umf)[3] <- time2sc

#recheck
summary( umf )
### end data prep -----------
### Analyze data ------------------------------------------
# We are now ready to perform our analysis. We start with a full model:
fm.closed <- pcount( ~ 1 + obsv + time + time2
                   ~ 1 + sagebrush + cheatgrass, 
                   #Define the maximum possible abundance
                   #during the primary occasion
                    K = 40,
                   mixture = c("P"),
                   data = umf )
# Note that we start with the observation submodel #
#We then define the ecological submodel 
# You should try alternative K values to make sure that your model 
#isn't sensitive to the value you provided. 

# View model results:
fm.closed

# Estimate confidence intervals:
confint( fm.closed, type = "state" )
# coefficients for detection submodel:
confint( fm.closed, type = 'det' )
#
# Based on the overlap of the 95% CIs for your predictor coefficients, #
# which may be important to each of your responses? #
# Answer:
#
#What is the mean abundance from our model?
exp(coef(fm.closed[1]))
plogis(coef(fm.closed[2]))
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
gof.boot <- Nmix.gof.test( fm.closed, nsim = 100, 
                           print.table = TRUE )
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
            K = 40, data = umf )
#view
fm.null
# Now build the list with the two models of interest:
rms <- fitList( 'full' = fm.closed,
                'null' = fm.null )
# Then use model selection function from unmarked, defining which is the null:
unmarked::modSel(rms, nullmod = "null" )
# What does this tell us?
# Answer:
#
# Now use gof checks outlined in Knape et al. 2018 MEE 9:2102-2114
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
# Niw plot Qq plots of randomized residuals against standard normal quantiles. #
# Under a good fit, residuals should be close to the identity line. 
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
# Are there other ways to fit time that may make more sense?
# Answer:
#

#########################################################################
##### Summarizing model output ##############
# Estimate partial prediction plots for predictors with 95% CIs not overlapping zero:
# Start by creating our datasets to predict over
# how many values do we use:
n <- 100
# what are the min max times:
closeddf %>% select( time.j1, time.j2, time.j3 ) %>% 
  summarise_all(list(min, max) )
#use them to define your bounds:
Time <- round(seq( 0, 360, length.out = n ),0)
Time2 <- Time^2
time.std <- scale( Time )
time2.std <- scale( Time2 )
# now for detection
detData <- data.frame( obsv = factor(c("tech.1", "tech.1","tech.1", "tech.1"), 
              levels = c("tech.1", "tech.2","tech.3", "tech.4") ), 
              time = time.std, time2 = time2.std )
#predict partial relationship:
pred.time <- predict( fm.closed, type = "det", newdata = detData, 
                           appendData = TRUE )

# create plot for detection submodel:
timep <- cbind( pred.time[,c("Predicted", "lower", "upper") ], Time ) %>%
  # define x and y values
  ggplot(., aes( x = Time, y = Predicted ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Time (mins pass 6:00am)", y = "Probability of Occupancy" ) +
  # add band of confidence intervals
  geom_smooth( aes(ymin = lower, ymax = upper ), 
               stat = "identity",
               size = 1.5, alpha = 0.5, color = "grey" ) +
  # add mean line on top
  geom_line( size = 2 ) 
#view
timep

############################################################################
################## Save your data and workspace ###################
# Save workspace:
save.image( "CountResults.RData" )

########## End of saving section ##################################

############# END OF SCRIPT #####################################