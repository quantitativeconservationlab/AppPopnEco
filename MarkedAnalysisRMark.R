#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we import our cleaned data for 2011 for single season capture-##
#  recapture analysis for Piute ground squirrels at the NCA.         ##
##                                                                   ##
# Squirrel abundance may be influenced by low temperatures when they ##
# emerge from hibernation in Feb and high temperatures in April-May, #
# and habitat composition including % of cheatgrass and sagebrush.    #
#                                                                     #
#                                                                     #            
# 20 sites were randomly selected for trapping over three days, after #
# three days of pre-baiting. This approach is meant to increase       #
# trappability, but may not avoid trap-happiness.                     #
# Trapping occurred over multiple years but tags used for marking     #
# individuals were temporary, so they lasted during the primary season#
# but not between seasons.                                            # 
#######################################################################

##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

# Start by installing the MARK program (which Rmark will talk to):
#MARK is freely available software that was written by Dr. Gary White  #
# and can be installed from:
# http://www.phidot.org/software/mark/index.html

#install RMark
install.packages( 'RMark' )

#Now load relevant packages:
library( tidyverse )
library( RMark ) 

## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Data/", sep = "" )
# load single season data:
closed_df <- read.csv( file = paste( datadir, "ind_2011.csv", sep = ""),
                    header = TRUE, colClasses = c("ch"="character") )
#view
head( closed_df ); dim( closed_df ) 
# load multi-season data:
open_df <- read.csv( file = paste( datadir, "ind_multi.csv", sep = ""),
                       header = TRUE, colClasses = c("ch"="character") )
#view
head( open_df ); dim( open_df ) 

#### End of data load -------------
####################################################################
##### Ready data for analysis --------------
# to use data with RMark we need  a column called ch that contains#
# the capture frequencies as a string, and a freq colum that #
# contains the capture frequencies (i.e., how many individuals with that
# observed capture history). If freq is absent (such as in our case, 
# then it assumed freq = 1).
# Grouping variables must be factors and individual covariates must #
# be numeric.

# let's check our dataframes
str(closed_df)
# We need to turn sex to factor 
closed_df$sex <- factor( closed_df$sex )

str(open_df)
# We need to turn sex to factor 
open_df$sex <- factor( open_df$sex )


huggins.df <- convert.inp( ind_df )
#sexid <- data.frame(sex=c("female","male")).

### end data prep -----------
### Analyze data ------------------------------------------
# We are now ready to perform our analysis. Model options available in RMark
# are in Table C.1 in the Laake, Rexstad 2008 Appendix C. 
# We concentrate on the following models: Closed, Huggins, Robust, RDHuggins
# Check out the manual for more details on each. 

# Single season models -----------
# The function mark is actually quite simple because it is a convenience function that calls 5 other
# functions that actually do the work in the following order:
#  1. process.data 
# 2. make.design.data 
# 3. make.mark.model 
# 4. run.mark.model
# 5. summary.mark
# We start with the single season models. Starting simply with a 
c.p <- process.data( closed_df, model = 'Closed', begin.time = 2011,
                  groups = "sex" )
names( c.p)
#define an intercept only model for p (or p[.] or M[0])
p.dot <- list( formula = ~1 )
# define p[sex]:
p.sex <- list(formula = ~sex )

# The parameter specifications are used with the mark argument model.parameters to define the
# model:
fm.c.0 <- mark( closed_df, model.parameters = )

####### comparing models
# If you leave this function empty, it searches at all objects with class 'mark'
# in the workspace and collates them into the ms object. 
ms <- collect.models()
# view 
ms
# if we need to remove some from the list 
#ms <- remove.mark( ms, c(1,3))
# if we wanted to model average results
ma <- model.average(ms,"p")
# Model averaging methods follow those in Burnham and Anderson (2002, Chpt 4)#
# CIs are estimated using the Delta-method. 
# Note that MARK can fail to adequately count the number of parameters in a #
# model, particularly complicated ones, and so it may favor overly complicated #
# models as a result of undercounting. So make sure that you check whether the #
# parameter count was done correctly. You can adjust by using adjust.parameter.count

# When would it be a good idea to model average?
# Answer:
#

# fixing parameters
p.time.fixed=list(formula=~time,fixed=list(time=c(2011,2018),value=0))

# You can set values to defaults if needed by:
#model0 <- mark( data, model.parameters = list(p=list(default = 0.9)))

# multi-season models ----------
# We start with the famous Cormack-Jolly-Seber (CJS) model #
# The model has two components: (1) a submodel for apparent survival (phi),
# and (2) a submodel for detection (p). See chpt 10 of Powell and Gale for #
# model details. 
# This model requires a single survey per Primary season, with #
# mortality allowed in between.
##########################################################################
# Model fit and evaluation -----------------------------------------------
# We start with goodness of fit (GoF) outlined by Duarte et al. 2018 #

#########################################################################
##### Summarizing model output ##############
#define the model that you want output for:
fm <- 
# The individual elements can be extracted using list notation. For example, 
#the data frame of the β parameters:
fm$results$beta
#beta estimates:
coefs(fm)
# To view all of the real parameters with standard errors, use summary
summary( fm, se = TRUE )
# Estimate partial prediction plots for predictors with 95% CIs not overlapping zero:
# Start by creating our datasets to predict over
# how many values do we use:
n <- 100
# Use the observed values to define our range:
cheatgrass <- seq( min( closeddf[,"cheatgrass"]),max( closeddf[,"cheatgrass"]),
                   length.out = n )
# what are the min max times:
closeddf %>% select( time.j1, time.j2, time.j3 ) %>% 
  summarise_all(list(min, max))
#use them to define your bounds:
Time <- round(seq( 0, 360, length.out = n ),0)
#standardize predictors:
cheat.std <- scale( cheatgrass )
time.std <- scale( Time )
#combine standardized predictors into a new dataframe to predict partial relationship
# for abundance submodel:
abundData <- data.frame( sagebrush = 0, cheatgrass = cheat.std )
#predict partial relationship:
pred.cheat <- predict( fm.closed, type = "state", newdata = abundData, 
                       appendData = TRUE )
#view
head( pred.cheat ); dim( pred.cheat )
# now for detection
detData <- data.frame( obsv = factor(c("tech.1", "tech.1","tech.1", "tech.1"), 
                                     levels = c("tech.1", "tech.2","tech.3", "tech.4") ), 
                       time = time.std )
#predict partial relationship:
pred.time <- predict( fm.closed, type = "det", newdata = detData, 
                      appendData = TRUE )

# create plot for ecological submodel:
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