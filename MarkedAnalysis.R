#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we import our cleaned data for single season capture-        ##
#  recapture analysis for Piute ground squirrels at the NCA.         ##
##                                                                   ##
## We start with analysis in unmarked. We will run 3 model options   ##
## M[0] = intercept only model for detection.                         # 
# M[t] = survey-varying detection model                              ##
# M[b] = behavior-effect on detection. To allow for trap-happiness or #
# trap-shyness following initial capture.                             #
#######################################################################

##### Set up your workspace and load relevant packages -----------
# install packages
install.packages( "AHMbook" )
#Now load relevant packages:
library( tidyverse )
library( unmarked ) 
library( AHMbook ) #contains data and functions from book

## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()
# load workspace with clean data
load( "MarkPrepWorkspace.RData")

#### End of data load -------------
####################################################################
##### Ready data for analysis --------------

# We need to custom build cell probability functions for each model type. 
# We start with M[t]
# Remember that custom crPiFun can only take a single argument p, which must be#
# the I x J matrix of detection probabilities. 
# we define the probabilities for each possible capture history:
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

# When providing a user defined piFun we also need information on how #
# to handle missing values by supplying a matrix of zeros (if not sampled) and ones (if sampled) #
# with rows = J (repeat surveys) and  cols = observed capture histories (ch):
o2y <- matrix( data = 1, nrow = J, ncol = CH )

# For the M[t] and M[0] models we can use our observations pooled #
# by capture histories. 
# A restriction of unmarked is that we cannot include #
# sites with all zero captures. i.e. sites where we trapped but #
# we caught no animals. We therefore work out which rows from our #
# dataframes we need to keep
keep <- which( rowSums(y_ik ) > 0 )
#define unmarked dataframe for M[t] and M[o] for single season
umf.Mt <- unmarkedFrameGMM( y = y_ik[keep,],
              #define abundance covariates
              siteCovs = as.data.frame(ik_sc[keep,]),
              #define observation covariates:
              obsCovs = list( wind = ij_wide[keep,widx],
                              temp = ij_wide[keep,tidx],
                              effort =  ij_wide[keep,eidx] ),
              #we use the pifun we defined for M[t]:
              obsToY = o2y, piFun = "crPiFun.t", 
              numPrimary = 1 )
#What does it look like?
summary(umf.Mt)

#Now define cell probabilities for the behavioral model M[b]:
#note this model assumes that on first capture, the animals are 'naive' #
# to the capture process. However, after being captured once, their #
# reaction to the capture process changes. 
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

### For the M[b] model we create a behavior dummy variable 
bMat <- matrix( c("Naive", "Wise", "Wise"), dim(y_ik[keep,])[1], 
                J, byrow = TRUE )
# Remember that there are limitations to the models that you can fit #
# for example, you cannot fit an M[b,t] model 

#define unmarked dataframe for M[b]
umf.Mb <- unmarkedFrameGMM( y = y_ik[keep,],
              #define abundance covariates
              siteCovs = as.data.frame(ik_sc[keep,]),
              #define observation covariates:
              obsCovs = list( behavior = bMat ),
              obsToY = o2y, piFun = "crPiFun.b", 
              numPrimary = 1 )
#check
summary(umf.Mb)
### end data prep -----------
#######################################################################
### Analyze data ------------------------------------------
# We are now ready to perform our analysis.
# For a single season: ######################
# We start with M[o] intercept only model to check it all runs: -------
M0.P <- gmultmix( #first function is abundance model. 
  #Second is availability. Third is detection.
            #we start with intercept only models
            ~1, ~1, ~1, 
          #we can provide the Mt dataframe
          umf.Mt, engine="R")
#Did it work?
M0.P
# We adjust to a Negative Binomial:
M0.NB <- gmultmix(   ~1, ~1, ~1, 
  umf.Mt, 
  mixture = "NB",
  engine="R" )
#Did it work?
M0.NB

#Now we run the full Mt model including also abundance predictors:
Mt.full.P  <- gmultmix( ~ 1 + shrub + annual + perennial,
                ~ 1 , 
                ~ 1 + wind + temp + effort,
                umf.Mt, 
                mixture = "P",
                engine="R" )

#Did it work?
Mt.full.P
#What do results say at first glance?
#Answer:
#
# We test whether the Negative Binomial as an alternative distribution
Mt.full.NB  <- gmultmix(  ~ 1 + shrub + annual + perennial,
                         ~ 1 , 
                         ~ 1 + wind + temp + effort,
                        umf.Mt, 
                        mixture = "NB",
                        engine="R" )

#check
Mt.full.NB

# We compare against the two full behavioral models
#The first uses the Poisson distribution:
Mb.full.P  <- gmultmix(  ~ 1 + shrub + annual + perennial,
                          ~ 1 , 
                          ~ 1 + behavior,
                          umf.Mb, 
                          mixture = "P",
                          engine="R" )
Mb.full.P
#The second full model uses the Negative Binomial:
Mb.full.NB  <- gmultmix(  ~ 1 + shrub + annual + perennial,
                          ~ 1 , 
                          ~ 1 + behavior,
                          umf.Mb, 
                          mixture = "NB",
                          engine="R" )

Mb.full.NB
# We therefore move on and try to see which of our Mt, Mo or Mb models #
# provided a better fit for our data: 
#Create model list:
Mtmodels <- fitList( "M0.NB" = M0.NB, 
                     "Mt.P" = Mt.full.P, 
                      "Mt.NB" = Mt.full.NB 
                       )

#compare
modSel( Mtmodels )
# Which is our top model?
# Answer:
# 
# What does it suggest: 
# Answer:
#
#Now repeat for Mb models
Mbmodels <- fitList( "M0.NB" = M0.NB, 
                     "Mb.P" = Mt.full.P, 
                     "Mb.NB" = Mt.full.NB 
)
modSel( Mbmodels )
##########################################################################
# Model fit and evaluation -----------------------------------------------
# We rely on unmarked options for boostrap goodness-of-fit option in #
# combination with custom built function fitstats, from Kery and Royle's
# Hierarchical modelling book, Chpt 7 pg 333. Also see example in ?parboot.
# The function computes error sum-of-squares, standard Chi-square and the #
# Freeman-Tukey statistic:
# Function returning three fit-statistics.
# fitstats <- function(fm) {
#   observed <- getY(fm@data)
#   expected <- fitted(fm)
#   resids <- residuals(fm)
#   sse <- sum(resids^2)
#   chisq <- sum((observed - expected)^2 / expected)
#   freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
#   out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
#   return(out)
# }
# the fitstats function is also available via the AHMbook package
# so we use that:
(gof.Mt.full.NB <- parboot( Mt.full.NB, fitstats2, nsim = 100) )

(gof.Mt.full.P <- parboot( Mt.full.P, fitstats2, nsim = 100) )

(gof.Mb.full.NB <- parboot( Mb.full.NB, fitstats2, nsim = 100) )

(gof.Mb.full.P <- parboot( Mb.full.P, fitstats2, nsim = 100) )

# What do these results suggest? 
# Answer: 
# 


#########################################################################
##### Summarizing model output ##############
# Select model to use
fm <- Mt.full.P
#Extract abundance estimates from top model
#start by creating a dataframe to hold results
outdf <- ik_df %>% select( idno, id, SiteID, year ) %>% 
        dplyr::filter( idno %in% rownames(y_ik[keep,]) )
#now calculate the mean and 95%CIs using the ranef function
outdf <- cbind( outdf, bup(ranef( fm)),confint( ranef( fm)))
colnames(outdf)[5:7] <- c( "Mean", "Lower", "Upper")
#Add raw counts
outdf$rawtrap <- as.vector(rowSums( y_ik[keep,] ))
#check
head(outdf)
#plot
ggplot( outdf, aes( x = as.factor(SiteID), y = Mean ) ) +
  theme_classic(base_size = 15) +
  geom_point( size = 3 ) +
  # add confidence intervals
  geom_errorbar( aes(ymin = Lower, ymax = Upper ), 
                 size = 1.5, width = 0.3 )  +
  facet_wrap( ~year, ncol = 1 )


# Estimate partial prediction plots for predictors with 95% CIs not overlapping zero:
# how many values do we use:
n <- 100
# Use the observed values to define our range:
shrub <- seq( min( ik_df[,"shrub"]), max( ik_df[,"shrub"]),
                   length.out = n )
#scale
shrub.std <- scale( shrub )
#combine standardized predictors into anew dataframe to predict partial relationship
# for abundance submodel:
abundData <- data.frame( shrub = shrub.std, annual = 0, 
                         perennial = 0  )
#predict partial relationship:
pred.shrub <- predict( fm, type = "lambda", newdata = abundData, 
                       appendData = TRUE )
#view
head( pred.shrub ); dim( pred.shrub )

# create plot for ecological submodel:
shrubp <- cbind( pred.shrub[,c("Predicted", "lower", "upper") ], shrub ) %>%
  # define x and y values
  ggplot(., aes( x = shrub, y = Predicted ) ) + 
  #choose preset look
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Shrub cover (%)", y = "Relative abundance" ) +
  # add band of confidence intervals
  geom_smooth( aes(ymin = lower, ymax = upper ), 
               stat = "identity",
               size = 1.5, alpha = 0.5, color = "grey" ) +
  # add mean line on top
  geom_line( size = 2 ) 
#view
shrubp
# How do you interpret this relationship?
# Answer:

#Compare abundance and partial prediction plot for the 
# model with the NB distribution instead. #
# What differences do you see? What distribution would you 
# choose moving forward? # 
# Answer: 

###################################################################
############################################################################
################## Save your data and workspace ###################
# Save workspace:
save.image( "MarkedResults.RData" )

########## End of saving section ##################################

############# END OF SCRIPT #####################################