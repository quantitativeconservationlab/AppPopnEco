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
library( unmarked ) #why unmarked?

## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Data/", sep = "" )
# load cleaned data
ind_df <- read.csv( file = paste( datadir, "ind_2011.csv", sep = ""),
                    header = TRUE )
#view
head( ind_df ); dim( ind_df ) 
#### End of data load -------------
####################################################################
##### Ready data for analysis --------------

huggins.df <- convert.inp( ind_df )
#sexid <- data.frame(sex=c("female","male")).

### end data prep -----------
### Analyze data ------------------------------------------
# We are now ready to perform our analysis. We start with a full model:

##########################################################################
# Model fit and evaluation -----------------------------------------------
# We start with goodness of fit (GoF) outlined by Duarte et al. 2018 #

#########################################################################
##### Summarizing model output ##############
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