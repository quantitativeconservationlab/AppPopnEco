install.packages( "R2ucare" )
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
##########################################################################

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
