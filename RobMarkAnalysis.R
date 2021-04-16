#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we import our cleaned data for 2011 for single season capture-##
#  recapture analysis for Piute ground squirrels at the NCA.         ##
##                                                                   ##
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

#install packages
install.packages( "R2ucare" )
#load packages
library( tidyverse )
library( RMark )
library( R2ucare )

## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Data/", sep = "" )
# load multi-season data:
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
# We need to turn sex to factor 
open_df$sex <- factor( open_df$sex )
#convert sites to factor
open_df$o.sites <- factor( open_df$o.sites)
#convert year to factor
open_df$year <- factor( open_df$year )

### Analyze data ------------------------------------------
# We are now ready to perform our analysis. Model options available in RMark
# are in Table C.1 in the Laake, Rexstad 2008 Appendix C. 
# We concentrate on the following models: Robust, RDHuggins
# Check out the manual for more details on each. 

# multi-season models ----------
# We start with the famous Cormack-Jolly-Seber (CJS) model #
# The model has two components: (1) a submodel for apparent survival (phi),
# and (2) a submodel for detection (p). See chpt 10 of Powell and Gale for #
# model details. 
# This model requires a single survey per Primary season, with #
# mortality allowed in between.


# The function mark is a convenience function that calls 5 other
# functions that do the following:
# 1. process.data 
# 2. make.design.data 
# 3. make.mark.model 
# 4. run.mark.model
# 5. summary.mark

# However, as you saw in Laake, Rexstad (2008) we can also manipuate each #
# step independently. 
# We start with step 1: process.data
# We need to specify the model type
o.pr <- process.data( open_df, model = 'Huggins',
                      #and how many groups are in our data
                      groups = c("sex", "o.sites", "year" ) )
summary( o.pr)
#Now we do step 2. make design data:
o.ddl <- make.design.data( o.pr ) 
# we have already incorporated predictors into our dataframe so we don't alter
# our ddl here
# Check that it worked as we expected:
o.ddl

# Now we define the formula of the submodels that we are interested in:
# Starting with the M[0] or M[.] or intercept-only model:
dot <- list( formula = ~1 )
# For detection we are interested in sex differences 
sex <- list( formula = ~ sex -1 )
sex.year <- list( formula = ~ sex  + year  -1 )
sex.site <- list( formula = ~ sex + o.sites  -1 )

fm.dot <- mark( o.pr, o.ddl, 
                model.parameters = list( p = dot, c = dot ))
fm.sex  <- mark( o.pr, o.ddl, 
                  model.parameters = list( p = sex, c = sex ))
# Now we run our competing model:
fm.sex.site  <- mark( o.pr, o.ddl, 
                  model.parameters = list( p = sex.site, c = sex.site ))
fm.sex.year <- mark( o.pr, o.ddl, 
                     model.parameters = list( p = sex.year, c = sex.year ))

# View results
fm.site$results$derived

##########################################################################
####### comparing models
# If you leave this function empty, it searches at all objects with class 'mark'
# in the workspace and collates them into the ms object. 
ms <- collect.models()
# view 
ms
# What do these results tell us?
# Answer:
# 

##########################################################################
# Model fit and evaluation -----------------------------------------------
# No goodness of fit options  for model types available in Rmark. A 
# median c-hat is available but only for limited models. Models cannot #
# have individidual covariates. See more details here:
# https://sites.warnercnr.colostate.edu/gwhite/median-chat/
# CJS models have goodness of fit options from package RELEASE:
# https://jamesepaterson.github.io/jamespatersonblog/2020-05-20_gof_for_CJS
#
# Model fit and evaluation -----------------------------------------------
# We start with goodness of fit from RELEASE:
# https://rdrr.io/cran/RMark/man/release.gof.html
# runs results for TEST2 and TEST3
# Test 2 = Does recapture depend on when an animal was first marked? Tests the equal catchability assumption.
# Test 3 = Does marking affect survival? Tests the equal survival assumption.
#https://jamesepaterson.github.io/jamespatersonblog/2020-05-20_gof_for_CJS

release.gof( c.pr )

#a matrix with a column for each capture event and a row for each individual
head( closed_df )
gofmat <- closed_df %>% select( trap.j1, trap.j2, trap.j3 )
gofmat <- as.matrix( gofmat)
colnames(gofmat) <- NULL
gofmat[1:10,]
table( closed_df$ch)
#Test 2.CT tests whether there is a difference in p at t+1 between those #
#captured and not captured at t when animals are known to be alive because #
# they are recaptured later in the study. In other words test for homegeneity #
# in captures
test2.t <- test2ct( gofmat[1:100,], rep(1, 100))# nrow(closed_df) ) )
# first argument = capture history matrix, second argument is frequency #
# of each capture history (1 for example) #

#Test2.CL tests if there is a difference in the expected time of next #
#recapture between individuals captured and not captured at t when animals #
# are known to be alive.
test2.l <- test2cl( gofmat, rep(1, nrow(gofmat) ) )

#Test 3 tests whether marking affects survival (equal survival assumption).#
#There are two components to Test 3 (Test3.SR and Test3.SM).
