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

#import predictor data
preddf <- read.csv( file = paste( datadir, "predictors.csv", sep = ""),
                    header = TRUE )
#view
head( preddf ); dim( preddf )
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

# let's check our dataframe
str(closed_df)
#remove unnecesary columns
c_df <- closed_df %>% select(o.sites, sex, ch )
#combine with site-level predictors:
c_df <- preddf %>% filter( marked == "yes" & year == 2011 ) %>% 
            #remove weather variables:
            select( o.sites, cheatgrass, sagebrush ) %>% 
            # right join so that only sites with captured individuals are kept:
            right_join( c_df, by = "o.sites" )
# We need to turn sex to factor 
c_df$sex <- factor( c_df$sex )
#convert sites to factor
c_df$o.sites <- factor( c_df$o.sites)
#view
head( c_df)
# now scale predictors
c_df$sagebrush <- scale( c_df$sagebrush )
c_df$cheatgrass <- scale( c_df$cheatgrass )
### end data prep -----------
### Analyze data ------------------------------------------
# We are now ready to perform our analysis. Model options available in RMark
# are in Table C.1 in the Laake, Rexstad 2008 Appendix C. 
# We concentrate on the following models: Huggins and Closed 
# Check out Lukacs Chapter 14 for more details on each. 
# What are the main differences between the two model options:
# Answer:
# 
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
c.pr <- process.data( c_df, model = 'Huggins',
                      #and how many groups are in our data
                      groups = c("sex", "o.sites" ) )
summary( c.pr)
#Now we do step 2. make design data:
c.ddl <- make.design.data( c.pr ) 
# we have already incorporated predictors into our dataframe so we don't alter
# our ddl here
# Check that it worked as we expected:
c.ddl

# Now we define the formula of the submodels that we are interested in:
# Starting with the M[0] or M[.] or intercept-only model:
dot <- list( formula = ~1 )
# For detection we are interested in sex differences and the effects of sagebrush:
sex.sage <- list( formula = ~sex + sagebrush -1 )
# Why do we have a -1?
# Answer:
# 
# What if differences in detection are due to site differences not related to 
# habitat? Could we add a site predictor instead? 
sex.site <- list( formula = ~sex + o.sites  -1 )
# Can you see that this is a competing model for the sex.sage option?

# Now we are ready for step 4: run our first model:
# The parameter specifications are used with the mark argument model.parameters 
#to define the model:
fm.sage <- mark( c.pr, c.ddl, 
          # the Huggins model has two response parameters: p = prob of capture
          #c = probability of recapture. We assign formulas to each submodel:
          model.parameters = list( p = sex.sage, c = sex.sage ))
#Step 5: summary of results:
summary( fm.sage)
#view coefficients:
fm.sage$results$beta
#view abundance estimates
fm.sage$results$derived
# What do these results tell us?
# Answer:
#
# Now we run our competing model:
fm.site  <- mark( c.pr, c.ddl, 
            model.parameters = list( p = sex.site, c = sex.site ))
# View results
fm.site$results$derived
# What do these tell us?
# Answer:

# Now let's run the null model to compare it against:
fm.null <- mark( c.pr, c.ddl, 
                 model.parameters = list( p = dot, c = dot ))
# view results:
fm.null$results$derived

#compare to naive estimates
colSums(closed_df[,c("trap.j1", "trap.j2", "trap.j3")])

####### comparing models
# If you leave this function empty, it searches at all objects with class 'mark'
# in the workspace and collates them into the ms object. 
ms <- collect.models()
# view 
ms
# What do these results tell us?
# Answer:
# 

# if we need to remove some from the list 
#ms <- remove.mark( ms, c(1,3))
# if we wanted to model average results
# ma <- model.average(ms,"p")
# Model averaging methods follow those in Burnham and Anderson (2002, Chpt 4)#
# CIs are estimated using the Delta-method. 
# Note that MARK can fail to adequately count the number of parameters in a #
# model, particularly complicated ones, and so it may favor overly complicated #
# models as a result of undercounting. So make sure that you check whether the #
# parameter count was done correctly. You can adjust by using adjust.parameter.count

# When would it be a good idea to model average?
# Answer:
#

# What about our 'Closed' model? ------------
c.pr.c <- process.data( c_df, model = 'Closed',
                      groups = c("sex", "o.sites" ) )
summary( c.pr.c)
#Now we do step 2. make design data:
c.ddl.c <- make.design.data( c.pr.c ) 
#step 4
fm.c <- mark( c.pr.c, c.ddl.c, 
                 model.parameters = list( p = sex.sage, c = sex.sage, f0 = dot ))
# View results
fm.c$results$derived


# Additional options
# fixing parameters
#p.time.fixed=list(formula=~time,fixed=list(time=c(2011,2018),value=0))

# You can set values to defaults if needed by:
#model0 <- mark( data, model.parameters = list(p=list(default = 0.9)))


#########################################################################
# Summarizing results
# Plot some partial prediction plots
# Answer:

############################################################################
################## Save your data and workspace ###################
# Save workspace:
save.image( "MarkResultsRMark.RData" )

######### End of saving section ##################################

############# END OF SCRIPT #####################################