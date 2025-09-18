#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we evaluate model fit and plot model output.                 ##
# After we conduct an analysis we get results for the single model.   #
# If we ran model selection, we get a top model. But is it any good?  #
# A model may over or under-fit the data and thus be a poor           #
# representation of the system we are trying to understand.           #
#  Too often ecologists stop when they finish running                 #
# their analyses. This is dangerous!                                  #
# Part of the issue is that evaluating models is hard. This is        #
# particularly so for Binomial models, as traditional approaches      #
#  that are suitable with Gaussian distributions do not work          #
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
# by estimating pseudo-R^2, based on Nagelkerke, N.J.D. (1991) A Note #
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
# For homework compare the full model against the 'top' model selected by 
#model selection last week here:
#

# Did the top model explain any more variation than the full model?
# Answer:
# 

# Note that there are multiple approaches for estimating pseudo-Rsquares. See:
# https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-pseudo-r-squareds/
# We chose the one available in unmarked but like everything, it should be #
# taken with a grain of salt. 


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
# For our occupancy submodel we have two continuous predictors so we can #
# create a partial plot of each while keeping the other predictor at the #
# mean value. 

# To create nice plots we first need to create new vectors with evenly #
# spaced values for our continuous predictors within their actual observed range: #
# Why do we not predict outside the observed range?
# Answer:
# 

# we start choosing the number of values we want to predict over:
n <- 100
# start with sagebrush
# we use our data (unscaled) to extract the observed range of the predictor:
sagebrush <- seq( min( closeddf[,"sagebrush"]),max( closeddf[,"sagebrush"]),
                  length.out = n )
#view
sagebrush
#standardize these values
sage.std <- scale( sagebrush )
#view
sage.std

#combine standardized predictor with other predictors in the model, kept at #
# their standardized mean:
sageData <- data.frame( sagebrush = sage.std, cheatgrass = 0 )
# Note that you have to label the  columns exactly the same as the names of your
# predictors in the model

#Use predict function to predict expected probability of occupancy across the 
# range of values you determine above for your predictor of interest, while 
# keeping other predictors at their mean:
pred.occ.sage <- predict( fm.closed, type = "state", newdata = sageData, 
                          appendData = TRUE )
# Note that you have to define which submodel you want to use for the prediction
# using the type = 'state' 
#view results
head( pred.occ.sage ); dim( pred.occ.sage )

### now replicate the process of predicting occupancy for observed values
# of cheatgrass here:
#Answer: 
#

### In our detection submodel we have a combination of continuous and #
# categorical predictors (where one level is the intercept) so we need # 
# to choose a level of categorical if we want to predict over the continuous #
# predictor

#sagebrush is also in the detection model so we don't have to re-create a new 
# vector to use it here. We use the same vector and combine it into a new dataframe
# with the other predictors in the detection model. Because they are categorical, 
# we need to chooose a reference category. Here we set it as tech 4:
sageDet <- data.frame( obsv = factor( "tech.4", levels = c("tech.1", "tech.2",
                  "tech.3", "tech.4") ), sagebrush = sage.std )

#view
head( sageDet )
#can you see that we still need to provide all levels because it is a factor #
# even though only one is chosen as the output?

#predict partial relationship between sagebrush and detection probability 
# using the predict function
pred.det.sage <- predict( fm.closed, type = "det", newdata = sageDet, 
                          appendData = TRUE )
# note that we define which submodel we want to use with type = 'det

#now if we want to look at differences in detection for the different levels
# of the categorical variable, we create a dataframe that varies those levels
# and keep the other predictors in that submodel at their mean value:
obsvDet <- data.frame( obsv = factor( c("tech.1", "tech.2","tech.3", "tech.4"), 
                       levels = c("tech.1", "tech.2","tech.3", "tech.4") ), 
                      sagebrush = 0 )
#view
obsvDet
#Now predict partial relationship between observer effects and detection
pred.det.obsv <- predict( fm.closed, type = "det", newdata = obsvDet, 
                          appendData = TRUE )
#view
pred.det.obsv

####### We are ready to visualize our results  #######
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

#Now add the plot for cheatgrass:
# Answer:
#

# visualize partial plots for detection
# Start with sagebrush:
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

# Now observer effects:
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