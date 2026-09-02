#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we import our cleaned data for one season of our occurrence   ##
#  observations for ground squirrels at the NCA and run a      ##
## closed population occupancy analysis. See Mackenzie et al. 2002   ##
## for details of the model. The occupancy model is hierarchical with #
# two components: (1) an ecological submodel linking occupancy to    ##
## environmental predictors at the site. (2) an observation submodel ##
## linking our detection probability to relevant predictors.         ##
##                                                                   ##
#######################################################################
##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

# Install new packages from "CRAN" repository. # 
install.packages( "unmarked" ) #package for estimating occupancy, N-mixtures, 
#and some multinomial approaches for capture data
install.packages( "MuMIn") # package for model selection and evaluation
# load packages:
library( tidyverse )#includes dplyr, tidyr and ggplot2
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( unmarked ) #main package for frequentist population models
library( MuMIn ) #model validation package
## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Data/", sep = "" )
# load our cleaned, single-season data
closeddf <- read.csv( file = paste( datadir, "closedf.csv", sep = ""),
                     header = TRUE )
#view
head( closeddf ); dim( closeddf ) 
#### End of data load -------------
####################################################################
##### Ready data for analysis --------------
# the unmarked function has several functions to make data import easy
# We need to define which predictors we will link to which responses 
# (i.e, detection or occupancy)
# Predictors for detection:
# We expect detection to be influenced by observer effects, but it could also #
# be affected by amount of cover obstructing visibility (so potentially a #
# negative relationship with sagebrush). #
# Predictors for occupancy:
# We expect occupancy to be influenced by vegetation (sagebrush and cheatgrass) #
# Why don't we include temperature in this model for one season?
# Answer: 
#

# For Homework modify the code below by including the temperature metrics into
# the occupancy model  and compare the results of the original full model 
# with vegetation only, vs a model that includes vegetation and temperature
# Do not include temperature in the model selection section!!! #

# Let's define the unmarked dataframe:
# Start by defining which columns represent the response (observed occurrences)
umf <- unmarkedFrameOccu( y = as.matrix( closeddf[ ,c("pres.j1", "pres.j2", "pres.j3")]),
      # Define predictors at the site level:
      siteCovs = closeddf[ ,c("sagebrush", "cheatgrass")],
    # Define predictors at the survey level as a list:
    obsCovs = list( obsv = closeddf[ ,c("observer.j1", "observer.j2", "observer.j3")] ) ) 
# View summary of unmarked dataframe:
summary( umf )

#now scale ecological predictors:
sc <- apply( siteCovs(umf), MARGIN = 2, FUN = scale )
# We replace the predictors in our unmarked dataframe with the scaled values:
siteCovs( umf ) <- sc
# Why do we scale predictors?
# Answer:
#

# View summary of scaled unmarked dataframe:
summary( umf )
# What does it tell us?
# Answer:
#

### end data prep -----------
### Analyze data ------------------------------------------
# We are now ready to perform our analysis. Since the number of predictors #
# is reasonable for the sample size, and there were no issues with #
# correlation, we focus on a single full, additive model:
fm.closed <- occu( #define detection submodel:
          ~ 1 + obsv + sagebrush 
          #define occupancy submodel second
          ~ 1 + sagebrush + cheatgrass, 
          #define dataframe to be used
          data = umf )
# Note that we start with the observation submodel, linking it to the intercept # 
# and observer effect, obsv. We then define the ecological submodel as related #
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
# Which predictors are important for each of your responses? #
# Answer:
# 

#############end full model ###########
###### Model selection and  collinearity assessments --------------------
# Indiscriminate model selection has become popular in recent years. #
# We do not believe model dredging is a suitable approach. #
# There are instances when model selection is required,
# When is model selection suitable?
# Answer:
#

# We demonstrate how to run alternative models
### Do not include weather metrics

# One reason for running alternative models is high collinearity among
# predictors. 
# As explained in Valle et a. 2024, collinearity can occur even when 
# correlation is low.

# Model outputs from unmarked do not work well with collinearity metrics
# BUT we can still evaluate collinearity among model predictors by running 
# a simple, single level (flat) model

#fit model with continuous predictors, choosing just detections from one day only
modcheck <- glm( pres.j2 ~ sagebrush + cheatgrass,
               data = closeddf, family = binomial)
#calculate VIF
sort( car::vif(modcheck), decreasing = T ) 

# What values of collinearity did you get?
# Answer:
#
# Which values would have been a concern? 
# Answer:
#

# Check that the choice of response (day of detections) in this model
# does not influence your results
# Answer:
# 

# Model selection 
# Let's assume that your vegetation predictors were found to be collinear and 
# you decide to run alternative models to see which may better explain occupancy

# We start by manually running alternative models:
#include a null model
( fm.1 <- occu( ~ 1 ~ 1, data = umf ) )
( fm.2 <- occu( ~ 1 + obsv + sagebrush  ~ 1 + sagebrush, data = umf ) )
( fm.3 <- occu( ~ 1 + obsv + sagebrush ~ 1 + cheatgrass, data = umf ) )

# why do we leave the detection model untouched?
# Answer:
#

# Use unmarked function we create a list of model options:
fms <- fitList( 'psi(.)p(.)' = fm.1,
  'psi(sagebrush)p(obsv+sagebrush)' = fm.2,
                'psi(cheatgrass)p(obsv+sagebrush)' = fm.3 )
#Note this uses the traditional (.) format to signify an intercept only model.
# We use unmarked function modSel() to compare models using AIC:
unmarked::modSel( fms )

# Which model(s) was/were most supported? Justify your answer.
# Answer:
#

# What would our estimates of occupancy be if we had not done any modeling?
# calculate naive occupancy by assigning a site as occupied if occurrence was #
# detected in any of the surveys, and as empty if ocurrence was not detected #
# in any of the surveys:
y.naive <- ifelse( rowSums( closeddf[ ,c("pres.j1", "pres.j2", "pres.j3")])>0,1,0)

# What are the estimates of occupancy from our models:
# Calculate Best Unbiased Predictors of site occupancy from each model:
# Estimate conditional occupancy at each site:
re <- ranef( fm.closed )
# the use those to estimate occupancy with the bup() function:
y.est.fm.closed <-round( bup(re, stat="mean" ) ) # Posterior mean
# Repeat this process for other models
y.est.fm.2 <-round( bup(ranef(fm.2), stat="mean" ) ) # Posterior mean
y.est.fm.3 <-round( bup(ranef(fm.3), stat="mean" ) ) # Posterior mean
# Compare results among them:
y.est.fm.closed - y.naive
y.est.fm.closed - y.est.fm.2
y.est.fm.closed - y.est.fm.3

# What do these results tell us?
# Answer:

# What was the estimated mean occupancy while keeping #
# sagebrush and cheatgrass at their mean values:
backTransform( linearComb( fm.closed, coefficients = c(1,0,0) , 
                           type = "state" ) )
# Note we transform the occupancy response (defined as state by unmarked) back #
# from the logit scale. The ecological model has 1 intercept and two predictors.#
# The predictors are scaled so their mean is 0, the intercept is 1, thus: c(1,0,0).#
# What was our estimated occupancy?
# Answer:
#
# What about our mean probability of detection for each observer?
# We start with observer 1:
backTransform( linearComb( fm.closed, coefficients = c(1,0,0,0,0), type = "det" ) )
#observer 2:
backTransform( linearComb( fm.closed, coefficients = c(1,1,0,0,0), type = "det" ) )
#observer 3:
backTransform( linearComb( fm.closed, coefficients = c(1,0,1,0,0), type = "det" ) )
#observer 4:
backTransform( linearComb( fm.closed, coefficients = c(1,0,0,1,0), type = "det" ) )
#mean occupancy for obsv 1 at mean % sagebrush:
backTransform( linearComb( fm.closed, coefficients = c(1,0,0,0,1), type = "det" ) )

# What do these results tell us about what drives occupancy and detection of #
# ground squirrels in your year of sampling?
# Answer:
#

# end of analysis ######

############################################################################
################## Save your data and workspace ###################

# This time we want to save our workspace so that we have access to all #
# the objects that we created during our analyses. #
save.image( "OccAnalysisWorkspace.RData" )

# Why don't we want to re-run the analyses every time instead?
# Answer:
#

########## End of saving section ##################################

############# END OF SCRIPT #####################################