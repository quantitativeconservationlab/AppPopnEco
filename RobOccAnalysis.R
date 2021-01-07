#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we import cleaned data under a robust design for occurrence  ##
#  observations of Piute ground squirrels at the NCA. We run a      ##
## dynamic occupancy analysis. See Mackenzie et al. (2003) Ecology   ##
## for details of the model. The dynamic occupancy model is hierarchical #
# with two components: (1) an ecological submodel that describes how ##
## site colonization and extinction processes explain site occupancy  ##
## with the opportunity to link these to environmental predictors at ##
## the site. (2) describes the observation submodel linking detection ##            
## probability to relevant predictors.                                ##
##                                                                   ##
#######################################################################
##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

# Install new packages from "CRAN" repository. # 

# load packages:
library( tidyverse )#includes dplyr, tidyr and ggplot2
library( unmarked ) #
## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Data/", sep = "" )
# load cleaned data with robust format
robdf <- read.csv( file = paste( datadir, "opendf.csv", sep = ""),
                      header = TRUE )
#view
head( robdf ); dim( robdf ) 
#### End of data load -------------
####################################################################
##### Ready data for analysis ----------------------------
# What predictors do we think drive colonization, extinction and # 
# detection of Piute ground squirrels at the NCA? #
# We build on some of our findings from the single-season occupancy #
# model. We found that cheatgrass and sagebrush were both potentially #
# important to occupancy. Breaking those relationships down into the #
# components of our dynamic process we predict that:
# Extinction:
# Cheatgrass is an invasive species that increases the likelihood of # 
# more frequent fires in the system so it may increase the probability #
# of a site becoming extinct. We also assess effects of min temperatures # 
# in Feb, assuming that colder Feb temperatures lower survival of adults #
# and reproduction if it limits food availability when ground squirrels #
# come out of hibernation. Particularly sites with low abundances, may become #
# extinct as a result. 
# Persistence:
# Sites with more sagebrush are more likely to persist as sagebrush #
#  provides food, and cover from the elements, and from predation. #

# Detection:
# We expect observer effects influence detection. We also test the #
# effects of sagebrush again. #

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