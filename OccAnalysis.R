#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we import our cleaned data for season 1 of our occurrence    ##
#  observations for Piute ground squirrels at the NCA and run a      ##
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
install.packages( "unmarked" ) #actually a collection of packages 

# load packages:
library( tidyverse )
library( unmarked )
## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Data/", sep = "" )
# load our cleaned data
closeddf <- read.csv( file = paste( datadir, "closedf.csv", sep = ""),
                     header = TRUE )
#view
head( closeddf ); dim( closeddf ) 
#### End of data load -------------
####################################################################
##### Ready data for analysis --------------
# the unmarked function has several functions to make data inport #
# easy
# We need to define which predictors we will link to which responses #
# We expect detection to be influenced by observer effects, but it could also #
# be affected by amount of cover obstructing visibility (so potentially a #
# negative relationship with sagebrush). #
# We expect occupancy to be influenced by habitat (sagebrush and cheatgrass) #
# Why don't we include temperature in this model for one season?
# Answer: 
#
# Let's define our unmarked dataframe:
# Start by defining which columns represent the response (observed occurrences)
umf <- unmarkedFrameOccu( y = as.matrix( closeddf[ ,c("pres.j1", "pres.j2", "pres.j3")]),
      # Define predictors at the site level:
      siteCovs = closeddf[ ,c("sagebrush", "cheatgrass")],
    # Define predictors at the survey level as a list:
    obsCovs = list( obsv = closeddf[ ,c("observer.j1", "observer.j2", "observer.j3")] ) ) 
#now scale ecological predictors:
sc <- apply( siteCovs(umf), MARGIN = 2, FUN = scale )
# We replace the predictors in our unmarked dataframe with the scaled values:
siteCovs( umf ) <- sc
# Why do we scale predictors?
# Answer:
#
# View summary of unmarked dataframe:
summary( umf )
# What does it tell us?

### end data prep -----------
### Analyze data ------------------------------------------
# We are now ready to perform our analysis. Since the number of predictors #
# is reasonable for the sample size, and there were no issues with #
# correlation, we focus on a single full, additive model:
fm.closed <- occu( ~ 1 + obsv + sagebrush
                   ~ 1 + sagebrush + cheatgrass, data = umf )
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
# How do we interpret these initial results?
# Answer:

# Indiscriminate model selection has become popular in recent years. #
# Although we do not believe this is a suitable approach here, we #
# demonstrate how to run various reduced models (not exhaustive list):
fm.2 <- occu( ~ 1 + obsv + sagebrush  ~ 1 + sagebrush, data = umf )
fm.3 <- occu( ~ 1 + obsv + sagebrush ~ 1 + cheatgrass, data = umf )
fm.4 <- occu( ~ 1 + obsv + sagebrush ~ 1, data = umf )
fm.5 <- occu( ~ 1 + obsv ~ 1, data = umf )
fm.6 <- occu( ~ 1 + sagebrush ~ 1, data = umf )
fm.7 <- occu( ~ 1 ~ 1, data = umf )
# Use unmarked function to create a list of model options:
fms <- fitList( 'psi(sagebrush + cheatgrass)p(obsv+sagebrush)' = fm.closed,
                'psi(sagebrush)p(obsv+sagebrush)' = fm.2,
                'psi(cheatgrass)p(obsv+sagebrush)' = fm.3,
                'psi(.)p(obsv+sagebrush)' = fm.4,
                'psi(.)p(obsv)' = fm.5,
                'psi(.)p(sagebrush)' = fm.6,
                'psi(.)p(.)' = fm.7 )
#Note this uses the traditional (.) format to signify an intercept only model
# We now compare models:
modSel(fms)
# Which model(s) was/were the most supported? 
# Answer:
#
# When would model selection be suitable?
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
# Repeat this process for other top model and the null:
y.est.fm.3 <-round( bup(ranef(fm.3), stat="mean" ) ) # Posterior mean
y.est.fm.7 <-round( bup(ranef(fm.7), stat="mean" ) ) # Posterior mean
# Compare results among them:
y.est.fm.closed - y.est.fm.3
y.est.fm.closed - y.est.fm.7
#view together
data.frame( y.naive, y.est.fm.closed, y.est.fm.3, y.est.fm.7 )
# What do these results tell us about the importance of accounting for effects #
# that impact detection?
# Answer:

# What was the estimated mean occupancy while keeping #
# sagebrush and cheatgrass at their mean values:
backTransform( linearComb( fm.closed, coefficients = c(1,0,0) , 
                           type = "state" ) )
# Note we transform the occupancy response (defined as state by unmarked) back #
# from the logit scale. The ecological model has 1 intercept and two predictors. The predictors are #
# scaled so their mean is 0, the intercept is 1, thus: c(1,0,0).
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
#  Piute ground squirrels in 2007?
# Answer:
#

# end of analysis ######
# Model fit and evaluation -----------------------------------------------
# After we conduct an analysis we also want to know if the model is any good.
# A model may be over or under-fit and thus a poor representation of the system #
# we are trying to understand. Too often ecologists stop when they finish running #
# their analyses. This is very dangerous. 

# One option to check the adequacy of model fit is to use a parametric bootstrap. #
# The parboot() function in unmarked simulate datasets from a fitted model, #
# refit the model, and generate a sampling distribution for a user-specified #
# fit-statistic #

# As recommended in the manual for binary data, we use a χ2 statistic: #
# Define the function that calculates the chi-squared statistic: 
chisq <- function(fm) {
  umf <- getData(fm)
  y <- getY(umf)
  y[y>1] <- 1
  sr <- fm@sitesRemoved
  if(length(sr)>0)
    y <- y[-sr,,drop=FALSE]
  fv <- fitted(fm, na.rm=TRUE)
  y[is.na(fv)] <- NA
  sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE)
}
# run
pb <- parboot(fm.closed, statistic = chisq, nsim = 100, parallel = FALSE)
#view results
pb


# calculate residuals
residuals( fm.closed )
plot( x = getY( fm.closed@data ), y = fitted( fm.closed ) )
#sum of squared residuals from model fit
SSE( fm.closed )

############################################################################

# if you want to save your workspace, because you are still #
# working through it use the following command:
save.image( "OccAnalysisWorkspace.RData" )


##########    save relevant data and workspaces     ###########
#save closed dataframe in our data folder:

########## End of saving section ##################################
################## Save your data and workspace ###################

############# END OF SCRIPT #####################################