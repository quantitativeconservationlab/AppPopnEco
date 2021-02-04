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
# in Feb, assuming that years with colder Feb temperatures lower survival of  #
# adults and reproduction if it limits food availability when ground squirrels #
# come out of hibernation. Particularly sites with low abundances, may become #
# extinct as a result. 
# Colonization:
# Sites with more sagebrush are more likely to become colonized as sagebrush #
# provides food, and cover from weather and predators. Years with cooler summers#
# may increase survival of young, and therefore the amount of individuals, #
# spreading to 

# Detection:
# We expect observer effects influence detection. We also test the #
# effects of sagebrush again. #

#determine year range
yrrange <- unique( robdf$year )
#years sampled (primary seasons)
T <- length( yrrange )
#sites sampled
M <- length( unique( robdf$o.sites ) )
#max replicate surveys per season 
J <- 3
J*T*M
#collapse wide presences and observer ids into long format:
longdf <- robdf %>% pivot_longer( cols = contains("j"), 
                                         names_to = c(".value", "survey"),
                                         names_pattern = "(.+).(.+)" )

#check
head( longdf );dim(longdf )
#do we have the right number of rows?
#answer: 

# Arrange columns in the correct order for unmarked fomatMult() function:
longdf <- longdf %>% select( year, o.sites, survey, pres., 
                             cheatgrass, sagebrush, Feb.minT, AprMay.maxT,
                             observer. ) 
# Use formatMulti to get data in correct format for robust analysis:
datadf <- formatMult( longdf )
#view
head( datadf )

#extract observations:
yy <- robdf %>%  
  select( o.sites, year, pres.j1, pres.j2, pres.j3 )
#check
head( yy )
#yyy <- matrix( as.character(yrrange), nrow(yy), T, byrow = TRUE )
#convert to wide format where each column is a yearXsurvey combination:
# we use pivot_wider instead of spread() so that we can do multiple columns
# at the same time:
yyy <- yy %>% 
  pivot_wider( names_from = year, #which criteria is getting turned to wide
               values_from = c("pres.j1", "pres.j2", "pres.j3" ), #what columns do we include
               #this rearranges the naming of the columns to use year first 
               names_glue = "{year}_{.value}"
               )
#check did we get the right number of dimensions?
head( yyy ); dim( yyy )
#now we sort columns by year instead of survey
yyy <- yyy %>% dplyr::select( str_sort( colnames(yyy), numeric = TRUE ) ) %>%
  # we remove site id
      select( -o.sites )

# our habitat occurs at the site level so summarize it accordingly:
sitecovs <- robdf %>% 
  group_by( o.sites ) %>%
  summarise( sagebrush = first( sagebrush ),
             cheatgrass = first( cheatgrass) ) %>%
  select( sagebrush, cheatgrass )
#check if it contains the right number of dimensions
head( sitecovs ); dim( sitecovs )

# our temperature predictors are at the year level but unmarked only has yearXsite 
# so we extract them that way:
Feb.minT <- robdf %>% 
  dplyr::select(o.sites, year, Feb.minT ) %>% 
  group_by(o.sites, year ) %>%
  dplyr::summarise( Feb.minT = first( Feb.minT ) ) %>%
    spread( key = year, value = Feb.minT )%>% 
  ungroup() %>% 
  dplyr::select( -o.sites )

#turn to matrix
Feb.minT <-as.matrix( Feb.minT )

#now for April-May temperatures 
AprMay.maxT <- robdf %>% 
  dplyr::select(o.sites, year, AprMay.maxT ) %>% 
  group_by(o.sites, year ) %>%
  dplyr::summarise( AprMay.maxT = first( AprMay.maxT ) ) %>%
  spread( key = year, value = AprMay.maxT ) %>% 
  ungroup() %>% 
  dplyr::select( -o.sites )
#turn to matrix
AprMay.maxT <- as.matrix( AprMay.maxT )

# our observer effects are at the observer-level and so we need to convert them #
# as such
obsv <- robdf %>%
        dplyr::select( o.sites, year,observer.j1, observer.j2, observer.j3 ) %>%
        pivot_longer( cols = starts_with("obs" ),
                      names_to = "survey", 
                      values_to = "obsv",
                      values_drop_na = FALSE )
#check
head( obsv ); dim( obsv )

trial <- robdf %>% 
  pivot_longer( everything(),
                names_to = ".value",
                names_pattern = "(.)(.)" )

# Let's define our unmarked dataframe:
# Start by defining which columns represent the response (observed occurrences)
umf <- unmarkedMultFrame( y = as.matrix( yyy ),
          # Define predictors at the site level:
          siteCovs = sitecovs,
          #define predictors at the year level:
          yearlySiteCovs = robdf[ ,c("Feb.minT", "AprMay.maxT") ],
          #yearlySiteCovs = list( Feb.minT = Feb.minT, AprMay.maxT = AprMay.maxT ),
          # Define predictors at the survey level as a list:
          #obsCovs = list( obsv = obsv[ ,"obsv"] ),
          # define the number of years
          numPrimary = T )


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
save.image( "RobOccAnalysisWorkspace.RData" )

# Why don't we want to re-run the analyses every time instead?
# Answer:
#

########## End of saving section ##################################

############# END OF SCRIPT #####################################