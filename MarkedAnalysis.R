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
# install packages
install.packages( "AHMbook" )
#Now load relevant packages:
library( tidyverse )
library( unmarked ) 
library( AHMbook ) #contains data and functions from book

## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Data/", sep = "" )
#load multi-season data 
open_df <- read.csv( file = paste( datadir, "ind_multi.csv", sep = ""),
                       header = TRUE, colClasses = c("ch"="character") )
#view
head( open_df ); dim( open_df ) 
# since we can model abundance we want to keep all sites in the single
# season analysis even if there were no captures. The closed_df we #
# had prepared in the MarkedPrep.R script removed those (for RMark)

#import predictor data
preddf <- read.csv( file = paste( datadir, "predictors.csv", sep = ""),
                    header = TRUE )
#view
head( preddf ); dim( preddf )
#### End of data load -------------
####################################################################
##### Ready data for analysis --------------
# start with the predictor dataframe: 
# we start by reducing our dataframe to our marked sites:
preddf <- preddf %>% filter( marked == "yes" & year == 2011 ) %>% 
  select( -counted, -marked, -yearname )
# create scaled dataframe: 
sc <- apply( preddf[ ,c("cheatgrass", "sagebrush")],
             MARGIN = 2, FUN = scale )
#turn to dataframe for input into umarked:
sitepreds <- as.data.frame( sc )
# check
head( sitepreds);dim( sitepreds)

#### now for our response dataframes: ----------
#define parameters
#number of replicate surveys
J <- 3
#number of sites
M <- length( unique( open_df$o.sites) )
# We can include abundance in these models, so we don't use #
# the closed_df we created in our MarkedPrep.R because that one #
# removed the sites with no captures. Here we keep them when we #
# assign o.sites as a factor prior to filtering:
#convert sites to factor
open_df$o.sites <- factor( open_df$o.sites )
#turn ch to factor for multi-season:
open_df$ch  <- factor( open_df$ch, 
                       levels=c("001", "010", "011", "100", "101", "110", "111"))

#number of observed capture histories:
CH <- length(levels(open_df$ch))

#make closed dataframe
closed_df <- open_df %>% filter( year == 2011 )
#check
head( closed_df ); dim( closed_df )

# We cannot use individual covariates so we collapse capture
# histories for each site:
y.closed <- table( closed_df$o.sites, closed_df$ch )
#turn to matrix
class( y.closed ) <- "matrix"
dim(y.closed)
head(y.closed)
#we need a matrix of zeros (if not sampled) and ones (if sampled) #
# with rows = J (repeat surveys) and  cols = observed capture histories (ch):
o2y <- matrix( data = 1, nrow = J, ncol = CH )

#Define cell probabilities manually for M[t] model for CH = 7:
# Remember that custom crPiFun can only take a single argument p, which must be#
# the M x J matrix of detection probabilities. 
crPiFun.t <- function(p) {
  p1 <- p[,1] #detection for J survey 1
  p2 <- p[,2] #detection for J survey 2
  p3 <- p[,3] #detection for J survey 3
  cbind("001" = (1-p1) * (1-p2) * p3,
        "010" = (1-p1) * p2 * (1-p3),
        "011" = (1-p1) * p2 * p3,
        "100" = p1 * (1-p2) * (1-p3),
        "101" = p1 * (1-p2) * p3,
        "110" = p1 * p2 * (1-p3),
        "111" = p1 * p2 * p3)
}
# To allow for p to vary by J survey we also need a covariate that represents 
# each survey ID to be specified as an obsCov
jMat <- matrix( as.character(1:3), M, J, byrow = TRUE )

#Now define cell probabilities for M[b] for CH = 7:
crPiFun.b <- function(p) { 
  pNaive <- p[,1]
  pWise <- p[,2]
  cbind("001" = (1-pNaive) * (1-pNaive) * pNaive,
        "010" = (1-pNaive) * pNaive * (1-pWise),
        "011" = (1-pNaive) * pNaive * pWise,
        "100" = pNaive * (1-pWise) * (1-pWise),
        "101" = pNaive * (1-pWise) * pWise,
        "110" = pNaive * pWise * (1-pWise),
        "111" = pNaive * pWise * pWise)
}

# Remember that there are limitations to the models that you can fit #
# for example, you cannot fit an M[b,t] model 
# we create a behavior covariate same as we did for the Jmat
bMat <- matrix( c("Naive", "Wise", "Wise"), M, J, byrow = TRUE )

#define unmarked dataframe for M[t] and M[o] for single season
umf.2011.Mt <- unmarkedFrameMPois( y = y.closed,
              siteCovs = sitepreds,
              #define survey id:
              obsCovs = list( J = jMat ),
              #we use the pifun we defined for M[t]:
              obsToY = o2y, piFun = "crPiFun.t" )
#Why don't we include temperature predictors in the single season?
#Answer: 
#
#check
summary(umf.2011.Mt)
#define unmarked dataframe for M[b] and M[o] for single season
umf.2011.Mb <- unmarkedFrameMPois( y = y.closed,
                siteCovs = sitepreds,
                #define behavior covariate:
                obsCovs = list( behavior = bMat ),
                #we use the pifum we defined for M[b]
                obsToY = o2y, piFun = "crPiFun.b" )

### end data prep -----------
#######################################################################
### Analyze data ------------------------------------------
# We are now ready to perform our analysis.
# For a single season: ######################
# We start with M[o] intercept only model to check it all runs: -------
Mo.int.closed <- multinomPois( #first function is for detection, second for abundance
            #we start with intercept only models for both
          ~1 ~1, 
          #we can provide the Mt dataframe
          umf.2011.Mt, engine="R")
#Did it work?
Mo.int.closed

#Now we add site-level predictors to abundance submodel
Mo.full.closed  <- multinomPois(~1
                                #abundance
                                ~1 + cheatgrass + sagebrush, 
                                umf.2011.Mt, engine="R")
#Did it work?
Mo.full.closed
#What do results say at first glance?
#Answer:
#
# Now Mt with all abundance covariates:
( Mt.full.closed <- multinomPois( ~1 + J 
                                 ~1 + cheatgrass + sagebrush, 
                                 umf.2011.Mt, engine="R" ) ) 
# Mb with all covariates for abundance:
( Mb.full.closed <- multinomPois( ~ behavior #-1
                                  ~1 + cheatgrass + sagebrush, 
                                  umf.2011.Mb, engine="R" ) ) 

# Note that we cannot include individual covariates in unmarked

# We therefore move on and try to see which of our Mt, Mo or Mb models #
# provided a better fit for our data: 
#Create model list:
closedmodels <- fitList( "Mo" = Mo.full.closed, 
                         "Mt" = Mt.full.closed 
                        ,"Mb" =  Mb.full.closed 
                       )

#why can't we compare these models?
# Answer: 
#
#compare
modSel( closedmodels )
# Which is our top model?
# Answer:
# 
# What does it suggest: 
# Answer:
#

##########################################################################
# Model fit and evaluation -----------------------------------------------
# We rely on unmarked options for boostrap goodness-of-fit option in #
# combination with custom built function fitstats, from Kery and Royle's
# Hierarchical modelling book, Chpt 7 pg 333. Also see example in ?parboot.
# The function computes error sum-of-squares, standard Chi-square and the #
# Freeman-Tukey statistic:
# Function returning three fit-statistics.
fitstats <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids^2)
  chisq <- sum((observed - expected)^2 / expected)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}
# the fitstats function is also available via the AHMbook package
(gof.Mb.full.closed <- parboot( Mb.full.closed, fitstats2, nsim = 1000,
                                report = 1) )

# What do this results suggest? 
# Answer: 
# 

#The last column on that table is actually meant to be the Bayesian p value
# estimate so that values <0.1 and >0.9 are actually signs of bad fit. 
# This suggests that our model is doing poorly. What else could you do 
# to improve model fit?

# Answer: 
#

#########################################################################
##### Summarizing model output ##############
# Estimate partial prediction plots for predictors with 95% CIs not overlapping zero:
# Which model do we want to use:
fm <- Mb.full.closed
# how many values do we use:
n <- 100
# Use the observed values to define our range:
cheatgrass <- seq( min( preddf[,"cheatgrass"]),max( preddf[,"cheatgrass"]),
                   length.out = n )
#scale
cheat.std <- scale( cheatgrass )
#combine standardized predictors into anew dataframe to predict partial relationship
# for abundance submodel:
abundData <- data.frame( cheatgrass = cheat.std, sagebrush = 0  )
#predict partial relationship:
pred.cheat <- predict( fm, type = "state", newdata = abundData, 
                       appendData = TRUE )
#view
head( pred.cheat ); dim( pred.cheat )

# we use the observed values to define our range:
sage <- seq( min( preddf[,"sagebrush"]),max( preddf[,"sagebrush"]),
             length.out = n )
#standardize them
sage.std <- scale( sage )
#combine standardized predictor into a new dataframe to predict partial relationship
sageData <- data.frame( cheatgrass = 0, sagebrush = sage.std )
#predict partial relationship
pred.sage <- predict( fm, type = "state", newdata = sageData, 
                      appendData = TRUE )
#view
head( pred.sage ); dim( pred.sage )
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

# create plots for ecological submodels
sagep <- cbind( pred.sage[,c("Predicted", "lower", "upper") ], sage ) %>%
  # define x and y values
  ggplot(., aes( x = sage, y = Predicted ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Sagebrush (%)", y = "Relative abundance" ) +
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
# Now we plot detection for trap naive and wise individuals:
#now for observer effects:
obsvDet <- data.frame( behavior = factor( c("Naive", "Wise" ), 
                levels = c("Naive", "Wise") ) )
#predict partial relationship between observer effects and detection
pred.det.obsv <- predict( fm, type = "det", newdata = obsvDet, 
                          appendData = TRUE )

# Now observer and detection:
obsvp.det <- pred.det.obsv %>%
  # define x and y values
  ggplot(., aes( x = behavior, y = Predicted, color = behavior ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  #remove legend
  theme( legend.position = "none" ) +
  # add labels
  labs( x = "Behavior", y = "Probability of detection" ) +
  #add mean detection for each observer
  geom_point( size = 4 ) +
  # add confidence intervals
  geom_errorbar( aes(ymin = lower, ymax = upper ), 
                 size = 1.5, width = 0.3 ) 
#view
obsvp.det
#What do these results suggest? Are individuals trap-happy or shy?
# Answer:
#
########################################################################
#### Removal sampling example ----------------------------
# We rely heavily on example 1, section 7.5 of hierarchical book:
# go here for more details: https://sites.google.com/site/appliedhierarchicalmodeling/home
data( ovendata )
#view bit of data:
ovendata.list$data[11:20,]; dim( ovendata.list$data )  
#so there are 70 point count sites and 4 time removal periods, 3min each:
#now check out the covariates:
ovendata.list$covariates[1:10,]
#where UFC is under foliage cover and TRBA is basal area of large trees
#create time removal dataframe
ovenFrame <- unmarkedFrameMPois(
  # import time removal counts:
  y = ovendata.list$data,
  #import and scale site level covariates:
  siteCovs = as.data.frame(scale(ovendata.list$covariates[,-1])),
  #define pifun type: 
  type = "removal" )
summary(ovenFrame)
#fit full model:
fm.tr <- multinomPois( #define detection model:
  ~ ufc 
  #define abundance model, including an interaction term:
  ~ ufc + trba + ufc:trba, 
  #define data
  ovenFrame )
#view results
summary(fm.tr)
#### end of removal example ----------
###################################################################
############################################################################
################## Save your data and workspace ###################
# Save workspace:
save.image( "MarkedResults.RData" )

#save the plot objects you need for your presentation
#Cheatgrass x abundance plot:
tiff( 'Data/TrapBehavior.tiff',
      height = 10, width = 12, units = 'cm', 
      compression = "lzw", res = 400 )
#call the plot
obsvp.det
#end connection
dev.off()

########## End of saving section ##################################

############# END OF SCRIPT #####################################