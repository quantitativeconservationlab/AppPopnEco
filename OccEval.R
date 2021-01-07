#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we evaluate model fit and plot model output.                 ##
# After we conduct an analysis we get results for the single model.   #
# If we ran model selection, we get a top model. But is it any good?  #
# We cannot know until we evaluate it! A model may over or under-fit  #
# the data and thus be a poor representation of the system we are trying #
# to understand. Too often ecologists stop when they finish running   #
# their analyses. This is very dangerous.                             #
# Part of the issue is the difficulties in evaluating a model. This is #
# particularly the case for binomial models, as several traditional   #
# approaches that are suitable with Gaussian distributions do not work #
# with binomial responses.                                            #
# There is therefore no silver bullet. We recommend exploring multiple #
# metrics as a way of confirming that your model is fit for purpose.  #
##                                                                   ##
#######################################################################
##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

#and some multinomial approaches for capture data
install.packages( "AICcmodavg" ) #for assessing goodness of fit

# load packages:
library( tidyverse )#includes dplyr, tidyr and ggplot2
library( unmarked ) #
library( AICcmodavg )

## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Data/", sep = "" )
# load our analysis workspace containing all the objects we created #
load( "OccAnalysisWorkspace.RData" )

###########end of loading ######
##########################################################################
# Model fit and evaluation -----------------------------------------------

# We assess goodness of fit (GoF) on detection frequencies, which relies on a #
# Pearson chi-square to assess fit as suggested by Mackenzie and Bailey (2004) #
# J. Agr., Bio. & Env. Stats. 9: 300-318
# using AICmodavg package
gof.boot <- AICcmodavg::mb.gof.test( fm.closed, nsim = 1000 )
#view
gof.boot

# What does the table tell us about the comparison between observed and expected #
# frequency of ground squirrel site-level detection histories?
# Answer:
# 
# What is the c-hat value?
# Answer:
# 
# Note that values of c-hat > 1 indicate overdispersion (variance > mean), but #
# that values much higher than 1 (i.e., > 4) probably indicate lack-of-fit. #
# In cases of moderate overdispersion, one usually multiplies the #
# variance-covariance matrix of the estimates by c-hat inflating the SE’s of the#
# estimates (c-hat is also known as a variance inflation factor). #
# In cases of underdispersion (c-hat < 1), keep the value of c-hat to 1. #
# Note that values of c-hat << 1 can also indicate lack-of-fit. #

# Is our model over- or under-dispersed?
# Answer:
#
# We can also evaluate how well our full model did against the null model # 
# by estimating pseudo-R^2, based on Nagelkerke, N.J.D. (2004) A Note #
# on a General Definition of the Coefficient of Determination. Biometrika 78,#
# pp. 691-692.#
# (1 - R^2) represents the proportion of unexplained variation in the model
# We create a reduced model list with only our two models of interest:
rms <- fitList( 'psi(sagebrush + cheatgrass)p(obsv+sagebrush)' = fm.closed,
                'psi(.)p(.)' = fm.16 )
# Then use model selection function from unmarked but this time we define #
# which one our null model is:
unmarked::modSel(rms, nullmod = "psi(.)p(.)" )
# What does this tell us about the fit of our model?
# Answer:
#
# Note that there are multiple approaches for estimating pseudo-Rsquares. See:
# https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-pseudo-r-squareds/
# We chose the one available in unmarked but like everything, it should be #
# taken with a grain of salt. 

# We could also compare all the models we run against our null:
# modSel(fms, nullmod = "psi(.)p(.)" )


# We can also estimate additional fit statistics borrowing a function #
# from the unmarked manual under parboot() section:
#function
fitstats <- function( model.name ) {
  #estimate observed y
  observed <- getY(model.name@data)
  # estimate expected occupancy given the model
  expected <- fitted(model.name)
  #calculate residuals using non-parametric bootstrapping
  resids <- residuals(model.name,
                      method = "nonparboot")
  #sum of squared residuals from model fit
  sse <- sum(resids^2, na.rm = TRUE)
  #chi-squared statistic:
  chisq <- sum((observed - expected)^2 / expected,
               na.rm = TRUE)
  #free Tukeys test statistic:
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, 
                  na.rm = TRUE)
  #output from the function
  out <- c(SSE = sse,
           Chisq = chisq,
           freemanTukey = freeTuke)
  
  return(out)
  
}
# define the model that you want to calculate statistics for:
model.name <- fm.closed
# run 
pb <- unmarked::parboot( model.name,
              fitstats,
              nsim = 100,
              report = TRUE
             )
#view
pb
#plot distribution of estimates for each statistic:
par(mfrow = c(3,1))
plot(pb, xlab = c("SSE", "Chisq", "FT") )

######## end of model evaluation ##############
##### Producing model output ##############
# Now that we have evaluated the value of our model we can produce #
# some output. If our aim is inference on what drives occupancy and #
# detection of Piute ground squirrels at the NCA, we would want to #
# plot the relationships between our model predictors and responses #

# Using our global model we estimate partial prediction plots #
# for our predictors. These plot the relationship between the #
# response and predictor of interest, while keeping the remaining #
# predictors at their mean values. # 
# What is the mean of our predictors?
# Answer:
#
# To create nice plots we need to create new vectors with evenly #
# spaced predictors within their actual observed range: #
# how many values do we use:
n <- 100
# we use the observed values to define our range:
sagebrush <- seq( min( closeddf[,"sagebrush"]),max( closeddf[,"sagebrush"]),
                  length.out = n )
cheatgrass <- seq( min( closeddf[,"cheatgrass"]),max( closeddf[,"cheatgrass"]),
                   length.out = n )
#standardize them
sage.std <- scale( sagebrush )
cheat.std <- scale( cheatgrass )
#combine standardized predictor into a new dataframe to predict partial relationship
# with sagebrush. We replace value of other predictor with its mean
sageData <- data.frame( sagebrush = sage.std, cheatgrass = 0 )
#predict partial relationship between sagebrush and occupancy
pred.occ.sage <- predict( fm.closed, type = "state", newdata = sageData, 
                          appendData = TRUE )
#view
head( pred.occ.sage ); dim( pred.occ.sage )

#combine standardized predictor into a new dataframe to predict partial relationship
# with cheatgrass. We replace value of other predictor with its mean
cheatData <- data.frame( sagebrush = 0, cheatgrass = cheat.std )
#predict partial relationship between sagebrush and occupancy
pred.occ.cheat <- predict( fm.closed, type = "state", newdata = cheatData, 
                          appendData = TRUE )

#combine standardized predictor into a new dataframe to predict partial relationship
# with sagebrush. We set observer effect as tech 1
sageDet <- data.frame( obsv = factor( "tech.1", levels = c("tech.1", "tech.2",
                  "tech.3", "tech.4") ), sagebrush = sage.std )
#predict partial relationship between sagebrush and detection
pred.det.sage <- predict( fm.closed, type = "det", newdata = sageDet, 
                          appendData = TRUE )
#now for observer effects:
obsvDet <- data.frame( obsv = factor( c("tech.1", "tech.2","tech.3", "tech.4"), 
                       levels = c("tech.1", "tech.2","tech.3", "tech.4") ), 
                      sagebrush = 0 )
#predict partial relationship between observer effects and detection
pred.det.obsv <- predict( fm.closed, type = "det", newdata = obsvDet, 
                          appendData = TRUE )

# create plots
#Starting with sagebrush and occupancy:
# select the predicted values we want to plot and combine with unscaled predictor
sagep <- cbind( pred.occ.sage[,c("Predicted", "lower", "upper") ], sagebrush ) %>%
  # define x and y values
  ggplot(., aes( x = sagebrush, y = Predicted ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Sagebrush (%)", y = "Predicted occupancy" ) +
  # add band of confidence intervals
  geom_smooth( aes(ymin = lower, ymax = upper ), 
               stat = "identity",
              size = 1.5, alpha = 0.5, color = "grey" ) +
  # add mean line on top
  geom_line( size = 2 ) 
#view
sagep

#Now cheatgrass and occupancy:
# select the predicted values we want to plot and combine with unscaled predictor
cheatp <- cbind( pred.occ.cheat[,c("Predicted", "lower", "upper") ], cheatgrass ) %>%
  # define x and y values
  ggplot(., aes( x = cheatgrass, y = Predicted ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Cheatgrass (%)", y = "Predicted occupancy" ) +
  # add band of confidence intervals
  geom_smooth( aes(ymin = lower, ymax = upper ), 
               stat = "identity",
               size = 1.5, alpha = 0.5, color = "grey" ) +
  # add mean line on top
  geom_line( size = 2 ) 
#view
cheatp

# Now sagebrush and detection:
# select the predicted values we want to plot and combine with unscaled predictor
sagep.det <- cbind( pred.det.sage[,c("Predicted", "lower", "upper") ], sagebrush ) %>%
  # define x and y values
  ggplot(., aes( x = sagebrush, y = Predicted ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Sagebrush (%)", y = "Predicted detection" ) +
  # add band of confidence intervals
  geom_smooth( aes(ymin = lower, ymax = upper ), 
               stat = "identity",
               size = 1.5, alpha = 0.5, color = "grey" ) +
  # add mean line on top
  geom_line( size = 2 ) 
#view
sagep.det

# Now sagebrush and detection:
# select the predicted values we want to plot and combine with unscaled predictor
obsvp.det <- pred.det.obsv %>%
  # define x and y values
  ggplot(., aes( x = obsv, y = Predicted, color = obsv ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  #remove legend
  theme( legend.position = "none" ) +
  # add labels
  labs( x = "Observer", y = "Predicted detection" ) +
  #add mean detection for each observer
  geom_point( size = 4 ) +
  # add confidence intervals
  geom_errorbar( aes(ymin = lower, ymax = upper ), 
               size = 1.5, width = 0.3 ) 
#view
obsvp.det

# If you were reporting these in a manuscript, which (if any) would you #
# leave out? Why? 
# Answer:
#
#### end model results section #############
###################
############################################################################
################## Save your data and workspace ###################

#save the plot objects we created 
#start by calling the file where you will save it
tiff( 'Data/SageXOcc.tiff',
      height = 10, width = 12, units = 'cm', 
      compression = "lzw", res = 400 )
#call the plot
sagep
#end connection
dev.off()
#Now the cheatgrass x occupancy plot:
tiff( 'Data/CheatXOcc.tiff',
      height = 10, width = 12, units = 'cm', 
      compression = "lzw", res = 400 )
#call the plot
cheatp
#end connection
dev.off()
#save the sagebrush x detection plot
tiff( 'Data/SageXDet.tiff',
      height = 10, width = 12, units = 'cm', 
      compression = "lzw", res = 400 )
#call the plot
sagep.det
#end connection
dev.off()

#save the observer x detection plot
tiff( 'Data/ObsvXDet.tiff',
      height = 10, width = 10, units = 'cm', 
      compression = "lzw", res = 400 )
#call the plot
obsvp.det
#end connection
dev.off()

# We save the workspace with our new objects in case we want to extract/modify #
# them #
save.image( "OccResultsWorkspace.RData" )


########## End of saving section ##################################

############# END OF SCRIPT #####################################