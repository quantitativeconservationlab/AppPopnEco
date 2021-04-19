#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we import our cleaned data for multi season capture-recapture ##
#            analysis for Piute ground squirrels at the NCA.         ##
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
install.packages( "lmerTest" )
install.packages( "visreg" )
install.packages( "pbkrtest" )
install.packages( "MuMIn")
install.packages( "DHARMa")
#load packages
library( tidyverse )
library( RMark )
library( lmerTest ) #allows fitting of mixed effect models for better diagnostics
library( visreg ) #plotting of mixed effect models
library( pbkrtest ) #to estimate better p-values for mixed-effects models
library( MuMIn ) # for model evaluation of mixed-effects models 
library( DHARMa ) #diagnostics for residuals in mixed-effect models
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

### Analyze detection data in RMark ------------------------------------------
# We are now ready to perform our analysis. Model options available in RMark
# are in Table C.1 in the Laake, Rexstad 2008 Appendix C. 
# For multi-season data models include: CJS, Robust, RDHuggins
# For more details on multiseason options check out the MARK book at: 
# http://www.phidot.org/software/mark/docs/book/

# For our example, remember that we do not have permanent marks and so #
# marks only last for each season. We therefore rely on single season #
# models and include year and site as groups to account for potential #
# differences in detection among those. Note that ideally we would be #
# able to incorporate random effects since we are repeatedly sampling the #
# same sites each year. Further we would ideally estimate abundance and detection #
# submodels in the same model, as these parameters are correlated with each other.
# Since we cannot do that in MARK we start by modeling detection only and #
# use the mean estimates of abundance from the top model in a second step #
# that links them to our predictors of interest. #
# I would normally instead use a Bayesian framework to analyse these data. #


# Remember, from Laake, Rexstad (2008) that we can manipuate each RMark #
# step independently. 
# We start with step 1: process.data
# We need to specify the model type
o.pr <- process.data( open_df, model = 'Huggins',
                      # We specify sex, sites and year as groups 
                      groups = c("sex", "o.sites", "year" ) )
summary( o.pr)
#For step 2, we create the design data:
o.ddl <- make.design.data( o.pr ) 
# we have already incorporated predictors into our dataframe so we don't alter
# our ddl here
# Check that it worked as we expected:
o.ddl

# For step 3, we define the formulas:
# Starting with the M[.], the intercept-only model:
dot <- list( formula = ~1 )
# Sex differences only
sex <- list( formula = ~ sex -1 )
# Sex and year
sex.year <- list( formula = ~ sex  + year  -1 )
# Sex and site
sex.site <- list( formula = ~ sex + o.sites  -1 )
# Why don't we run a sex, site and year model?
# Answer:
# 
# Step 4, we run our model alternatives:
#M[.]:
fm.dot <- mark( o.pr, o.ddl, 
                model.parameters = list( p = dot, c = dot ))
#M[sex]:
fm.sex  <- mark( o.pr, o.ddl, 
                  model.parameters = list( p = sex, c = sex ))
#M[sex + site]:
fm.sex.site  <- mark( o.pr, o.ddl, 
                  model.parameters = list( p = sex.site, c = sex.site ))
#M[ sex + year]:
fm.sex.year <- mark( o.pr, o.ddl, 
                     model.parameters = list( p = sex.year, c = sex.year ))

# Step 5, view results
summary( fm.sex.year )
fm.sex$results$derived

####### comparing models
# If you leave this function empty, it searches at all objects with class 'mark'
# in the workspace and collates them into the ms object. 
ms <- collect.models()
# view 
ms
# What do these results tell us?
# Answer:
# 

# From these we will choose the M[sex] as our top model. Remember this model #
# also includes a behavioral effect of trapping (squirrels are trap-happy) #

##########################################################################
# Model fit and evaluation for detection models -----------------------------------------------
# Check out chapter 5 of the MARK book for more details on model evaluation #
# options for CJS and other multi-season models available in program mark.#
# You can do these in R using package R2ucare. See an example in:
# https://jamesepaterson.github.io/jamespatersonblog/2020-05-20_gof_for_CJS
# The package runs results for TEST2 and TEST3
# Test 2 = Does recapture depend on when an animal was first marked? #
# Tests the equal catchability assumption.
# Test 3 = Does marking affect survival? Tests the equal survival assumption.#
# There are other subcomponents eg:
#Test 2.CT checks for differences in p at j+1 between those that were #
#captured and not captured at j for animals known to be alive because #
# they were recaptured later in the study. In other words test for homogeneity #
#Test2.CL tests if there is a difference in the expected time of next #
#recapture between individuals captured and not captured at t when animals #
# are known to be alive.
#Test 3 tests whether marking affects survival (equal survival assumption).#
#There are two components to Test 3 (Test3.SR and Test3.SM).

# HOWEVER, these tests are not suitable for our single season scenario

# there is a way to calculate c hat for our top model:
fm.sex.c.hat <- fm.sex$results$deviance/fm.sex$results$deviance.df

# Besides that, we can also look at AIC between the top model and the null
# which at least tells us that the top model explains a lot more of the #
# variation. Without any other current options we move forward to the second #
# step of analysis, which links abundance estimates from the top model to #
# our predictor. #

##### end of model evaluation ######
#### summarizing results for detection models ###############################################
# 
# covariate.predictions( fm.sex, data = data.frame( sex = c("male", "female")) )
# What are some of the limitations of this two-step approach?
# Answer: 
#

# extract initial detection estimates 
p_df <- get.real( model = fm.sex, parameter = 'p', se = TRUE)
# select relevant columns 
p_df <- p_df %>%  select(o.sites,year,time,sex,estimate,lcl,ucl) 
#check
head(p_df, 20); dim( p_df )

# note that initial detection is the same for each j survey, year and site. 
# So we only plot differences between sexes:
p_df %>% filter( time == 1 & year == 2007 & o.sites == 8 ) %>% 
        ggplot(., aes( x = sex, y = estimate, color = sex ) ) +
  #choose preset look
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Sex", y = "Initial detection" ) +
  #add mean detection for each observer
  geom_point( size = 4 ) +
  # add confidence intervals
  geom_errorbar( aes(ymin = lcl, ymax = ucl ), 
                 size = 1.5, width = 0.3 )
        
head(p_df)

# Extract group labels from p_df:
N_df <- p_df %>% filter( time == 1 ) %>% 
        select( o.sites,year,sex )
#extract abundance estimates from model results object:
N_df <- cbind( N_df, fm.sex$results$derived )
#check
head( N_df ); dim( N_df )
  
#Combine estimates across sexes:
N_df <- N_df %>% group_by( o.sites, year ) %>% 
  #sum abundance for both sexes
  summarise( totalN = sum( `N Population Size.estimate` ) ) %>% 
  #turn year back to numeric
  mutate( year = as.numeric(as.character(year)),
          o.sites = as.numeric(as.character(o.sites))) %>% 
  ungroup() 

# Plot annual changes in abundance by site:
N_df %>%
  ggplot(., aes( x = year, y = totalN ) ) +
  theme_bw( base_size = 15 ) +
  # add labels
  labs( x = "Year", y = "Abundance estimate" ) +
  #these are mean abundance estimates:
  geom_point( size = 4 ) +
  geom_line( size = 2 ) +
  #plot each site separately 
  facet_wrap( ~o.sites, ncol = 4 )
####### end of detection models ##################################
######################################################################
####### Abundance analysis, model evaluation and result summaries #########
# start by linking estimated abundance with predictor data:
head( preddf ); head( N_df )
# Combine with predictor dataframe:
N_df <- preddf %>% select( -counted, -marked, -yearname ) %>% 
        right_join( N_df, by = c( "o.sites", "year" ) )
#check
head( N_df )
# scale predictors:
N_df[,c("cheatgrass", "sagebrush", "Feb.minT", "AprMay.maxT")] <- apply( 
  N_df[ ,c("cheatgrass", "sagebrush", "Feb.minT", "AprMay.maxT")],
               MARGIN = 2, FUN = scale )

# We start by running a full model including all fixed effects of interest #
# as well as a random intercept for sites. 
am1 <- lmerTest::lmer( totalN ~ cheatgrass + sagebrush + Feb.minT + AprMay.maxT +
                 ( 1|o.sites ),
               data = N_df )
#view results
summary( am1 )

#What does the random intercept account for?
# Answer:
# 

# To define which predictors are important we start by calculating 'better' #
# p values using the Kenward-Roger approximation. Traditional estimates of p #
# values are not valid with random effect models.
# So instead of the Zvalue that comes with glmmTMB or lme4 we rely on lmerTest
# and use model output with the anova function to extract 'Kenward-Roger p values'
anova( am1, type = 2, ddf = "Kenward-Roger")

# Validate results by estimating 95% CIs
# This first option uses the Wald method
confint( am1 ) 
# the second option provides more robust estimates through bootstrapping:
am1.cis <- confint( am1, method = 'boot' )# 'profile' )
round( am1.cis, 2 )
# How much do the two types of CIs differ?
# Answer:
# 

# We can also calculate marginal R^2 for random effects only (R2m) and the #
#  whole model including fixed and random (R2c) using the MuMIn package:
MuMIn::r.squaredGLMM( am1 )
# How do we interpret these results?
# Answer: 
# 
# Tip: think about how much the R2c and R2m differ, and how much of that #
# is due to the fixed effects.

library( DHARMa )
model_simres <- simulateResiduals( am1 )
#you can plot the results:
plot( model_simres )
# What do these diagnostics tell us?
# Answer:
#

#### view partial relationships with important predictors:
visreg( am1, "cheatgrass" )
visreg( am1, "Feb.minT" )

# What would be our conclusions based on these results?
# Answer:
#
# How are our estimates likely to be affected by using this 2-step approach?
# Answer:
#
###### end abundance models ###########################
################## Save your data and workspace ###################
# Save workspace:
save.image( "RobMarkResults.RData" )

######### End of saving section ##################################

############# END OF SCRIPT #####################################