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

#Now load relevant packages:
library( tidyverse )
library( unmarked ) 

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
preddf <- preddf %>% filter( marked == "yes" ) %>% 
  select( -counted, -marked, -yearname )
#we scale them
sc <- apply( preddf[ ,c("cheatgrass", "sagebrush", "Feb.minT", "AprMay.maxT")],
             MARGIN = 2, FUN = scale )
# create scaled dataframe: 
preddf.sc <- preddf
preddf.sc[ ,c("cheatgrass", "sagebrush", "Feb.minT", "AprMay.maxT")] <- sc
#check
head( preddf.sc)
#extract site covariates for 2011 and removing unmarked sites:
sitepreds <- preddf.sc %>% filter( year == 2011 ) %>% 
  select( cheatgrass, sagebrush )
head( sitepreds)
dim( sitepreds)

#now create siteXyear dataframes:
#for sagebrush:
sagebrush <- preddf.sc %>% select( o.sites, year, sagebrush ) %>% 
  spread( key = year, value = sagebrush ) %>% 
  select( -o.sites )
#check
head(sagebrush); dim(sagebrush)
#for cheatgrass:
cheatgrass <- preddf.sc %>% select( o.sites, year, cheatgrass ) %>% 
  spread( key = year, value = cheatgrass ) %>% 
  select( -o.sites )
#check
head(sagebrush); dim(sagebrush)

#for Feb temperature:
Feb.minT <- preddf.sc %>% select( o.sites, year, Feb.minT ) %>% 
  spread( key = year, value = Feb.minT ) %>% 
  select( -o.sites )
#check
head(Feb.minT); dim(Feb.minT)
#for Apr-May temperature:
AprMay.maxT <- preddf.sc %>% select( o.sites, year, AprMay.maxT ) %>% 
  spread( key = year, value = AprMay.maxT ) %>% 
  select( -o.sites )
#check
head(AprMay.maxT); dim(AprMay.maxT)

#### now for our response dataframes: ----------
#define parameters
#number of replicate surveys
J <- 3
#number of sites
M <- length( unique( open_df$o.sites) )
#year range
yrrange <- sort( unique( open_df$year))
#number of years
T <- length( yrrange)
#turn ch to factor for multi-season:
open_df$ch  <- factor( open_df$ch, 
                       levels=c("001", "010", "011", "100", "101", "110", "111"))
#convert sites to factor
open_df$o.sites <- factor( open_df$o.sites )

#number of observed capture histories:
CH <- length(levels(open_df$ch))
#check
str(open_df)
#how many observations (frequency) for each capture history?
table( open_df$ch)

#we also need a matrix of zeros and ones with rows = J (repeat surveys) and 
# cols = number of observed capture histories (ch):
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
jMat <- matrix( 1:3, M, J, byrow = TRUE )

#Now define cell probabilities for M[b] for CH = 7:
crPiFun.b <- function(p) { 
  pNaive <- p[,1]
  pWise <- p[,3]
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
bMat <- matrix( c("Naive", "Naive", "Wise"), M, J, byrow = TRUE )

# since we can include abundance in our models then we don't use #
# the closed_df we created in our MarkedPrep.R because that one #
# removed the sites where we had no captures. Instead here we keep #
# them by having assigned o.sites as a factor prior to filtering:
closed_df <- open_df %>% filter( year == 2011 )
#check
head( closed_df ); dim( closed_df )
# if we are not using individual covariates then we can collapse our capture
# histories for each site as follows
y.closed <- table( closed_df$o.sites, closed_df$ch )
#turn to matrix
class( y.closed ) <- "matrix"
dim(y.closed)
#define unmarked dataframe for M[t] and M[o] for single season
umf.2011.Mt <- unmarkedFrameMPois( y = y.closed,
              siteCovs = sitepreds,
              #define survey id:
              obsCovs = list( J = jMat ),
              obsToY = o2y, piFun = "crPiFun.t" )
#Why don't we include temperature predictors in the single season?
#Answer: 
#
#check
umf.2011.Mt
#define unmarked dataframe for M[b] and M[o] for single season
umf.2011.Mb <- unmarkedFrameMPois( y = y.closed,
                siteCovs = sitepreds,
                #define behavior covariate:
                obsCovs = list( behavior = bMat ),
                obsToY = o2y, piFun = "crPiFun.b" )

####### Multi-season data prep ---------------
# to extract data for multi-seasons we start by selecting the #
# capture histories
head(open_df )
# we calculate the frequencies for each site and year:
y.df <- table( open_df$o.sites, open_df$ch,open_df$year )
# then we convert from list to wide format
y.open <- y.df[,,1]
for( t in 2:T){
  y.open <- cbind( y.open, y.df[,,t] )
}
#check it
head( y.open ); dim( y.open)
# convert to matrix
class( y.open ) <- "matrix"
#now expand the o2y:
o2yGMM <- kronecker(diag(T), o2y)

#updated dataframe for behavior repeated over T years:
bMat.open <- replicate( T, bMat, simplify = FALSE)
bMat.open <- do.call(cbind.data.frame, bMat.open)
head( bMat.open); dim( bMat.open)
#now turn it into unmarked dataframe:
umf.open <- unmarkedFrameGMM( y = y.open, 
            # we don't have any site level covariates
            # We define our siteXyear covariates:
            yearlySiteCovs = list(  sagebrush = sagebrush,
            cheatgrass = cheatgrass,Feb.minT =Feb.minT,
            AprMay.maxT =  AprMay.maxT ),
            #define behavior covariates:
            obsCovs = list( behavior = bMat.open ),
            #user defined cell probabilities and years
            obsToY = o2yGMM, piFun = "crPiFun.b", numPrimary = T )
          
#view
summary( umf.open )
### end data prep -----------
### Analyze data ------------------------------------------
# We are now ready to perform our analysis.
# For a single season: ######################
# We start with M[o] models: -------
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
( Mb.full.closed <- multinomPois( ~ behavior -1
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
######## For multiple seasons -----------

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

(gof.Mt.full.closed <- parboot( Mt.full.closed, fitstats, nsim = 1000,
                                report = 1) )

# What do this results suggest? 
# Answer: 
# 

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
save.image( "MarkedResults.RData" )

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