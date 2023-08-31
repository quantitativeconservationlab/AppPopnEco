#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                  ###
##  Our study species is the Piute ground squirrel, Spermophilus     ## 
## mollis. Ground squirrels are widely distributed in sagebrush      ##
## habitats of the Great Basin and Columbia Plateau. We focus        ##
##  on their populations inside the birds of prey NCA.               ##
## In this script we prepare the data that we will use              ###
## to estimate site occupancy of Piute ground squirrels.             ## 
##  Site locations were randomly selected within the study area     ###
## ensuring that they were at least 800 m apart from each other     ##
## Sampling occurred during 2007-2018, with 3 repeat surveys each year. #
## Field surveys consisted of active searches for individuals, active burrows  #
## or fresh scats within a 200m radius of the site center point.     ##
##  Four technicians worked on the project. The sites visited each   ##
### year were randomly assigned to an observer.                      ##
##  A subsample of sites were also surveyed with walking transects   ##
##  and live-trapping, which we will use for later analyses of      ##
## abundance. Time of day was noted during the count surveys.        ##
##                                                                   ##
## For our occupancy data we create two data summaries:            ##
## 1) for one year to demonstrate use of closed population models.    #
## 2) for all years to asssess temporal changes using aROBUST DESIGN. #              ##
##                                                                   ## 
## Predictors:                                                       ##
# We use landcover data from the National Geospatial Data Asset (NGDA).#
# https://www.mrlc.gov. We downloaded sagebrush and annual           ##
# grasses (which we assumed are mainly cheatgrass) from 2007-2018    ##
# Data are available as a raster of 30 x 30 m cell resolution.       ##
# We summarized annual % cover surrounding our sites with a 200 m buffer #
# Habitat data are at the site X year level resolution.               #
# Climate data were downloaded from the Twins Falls Station           #
# including daily measures of min, mean and max Temperatures          # 
# and total precipitation. We extracted minimum Temperatures for      #
# Feb, when squirrels come out of hibernation and max temperatures    #
# for Apr-May, when consistently hot days may prevent them from       #
# foraging.                                                           #
#######################################################################
##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

# Install new packages from "CRAN" repository. # 
install.packages( "tidyverse" ) #actually a collection of packages 

# load packages:
library( tidyverse ) 
## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are. 
# I have it in a Data folder in my Rstudio project:
datadir <- paste( getwd(), "/Data/", sep = "" )

# load observed occurrences:
obs_df <- read.csv( file = paste( datadir, "obs_df.csv", sep = ""),
                     header = TRUE )
#view
head( obs_df ); dim(obs_df)
#load predictor data
preddf <- read.csv( file = paste( datadir, "predictors.csv", sep = ""),
                   header = TRUE )
#view
head( preddf ); dim( preddf )
########## end of data load ################################
#######################################################################
######## explore data #############
# Start viewing our observations dataframe:
head( obs_df ); dim( obs_df )
# Our observations dataframe metadata:
# o.sites = site id (I)
# counted = was the site also surveyed using point counts
# marked = was the site also trapped
# year = year of sampling (primary seasons) T
# pres.j1, pres.j2, pres.j3: observation of occurrence for j surveys 1:3
# count.j1, count.j2, count.j3 = point counts for j surveys 1:3
# observer.j1, observer.j2, observer.j3 = observer IDs for each j survey, site
# time.j1, time.j2, time.j3 = time of day that the point count were conducted #
# where 0 is 0 minutes pass 6am to 360 minutes at noon.

# Let's define some parameters
# How many sites were sampled?
I <- length( unique( obs_df$o.sites ) )
# What years were sampled (primary seasons)?
yrrange <- sort( unique( obs_df$year) )
#How many years (i.e., primary seasons )?
T <- length( yrrange )
# How many replicate surveys
J <- 3
#view
I; yrrange; T; J
# For our 1st set of data we create a closed population dataframe #
# that only includes the last year data, and columns that are #
# relevant. Note we don't include time of day because it is only #
# relevant to point counts.
    # Select obs_df:
closeddf <- obs_df %>% 
    #filter only rows for last year:
    dplyr::filter( year == yrrange[T] ) %>%
    #select desired columns to keep:
    dplyr::select( o.sites, year, pres.j1, pres.j2, pres.j3,
            observer.j1, observer.j2, observer.j3 ) 
#view resulting dataframe
head( closeddf ); dim( closeddf )

# How many detections each survey across our 100 sites?
colSums( closeddf[, c("pres.j1", "pres.j2","pres.j3")])

# We also check for missing values in the response #
colSums( is.na( closeddf[, c("pres.j1", "pres.j2","pres.j3")]) )
#none are present

#### Checking our predictors ------------
# What about our predictors?
# view
head( preddf ); dim( preddf )
# Before thinking of whether we can include predictors in our#
# models we should check their distribution and correlation #
# Why?
# We start by checking for outliers, skewed distribution etc #
# create a vector with predictor names
prednames <- c("cheatgrass", "sagebrush", "Feb.minT", "AprMay.maxT" )
# loop over each to create histograms for each predictor:
for( p in 1:length(prednames) ){
  # create an object with the ggplot so that you can display it 
  # in a loop 
  a <- ggplot( preddf ) + #choose your data
    theme_bw( base_size = 15 ) + #choose a preset theme
    labs( x = prednames[p] ) + #label x axis using our predictor names
    geom_histogram( aes( get(prednames[p]) ), bins= 10 ) #plot histogram
  # display your plot object
  print( a )
}
# What do you note? Any apparent issues with these data?
# Answer:
#
# Let's plot how predictors vary annually:
for( p in 1:length(prednames) ){
  # We can also incorporate site variability for habitat:
  bp<- ggplot( preddf ) +
    theme_bw( base_size = 15 ) + #choose a preset theme
    theme( legend.position = "none" ) + #remove legend display
    labs( y = prednames[p] ) + #label x axis using our predictor names
    # plot each site individually
    geom_line( aes( x = year, y = get(prednames[p]),
    color = as.factor(o.sites) ), size = 1.5 )
  
  # Here we rely on a smoothing spline to get mean annual trends #
  # across all sites:
  cp<- ggplot( preddf, aes( x = year, y = get(prednames[p]) ) ) +
    theme_bw( base_size = 15 ) + #choose a preset theme
    labs( y = prednames[p] ) + #label x axis using our predictor names
    geom_smooth( size = 2 ) #smooth mean across all sites
  #display plots
  print( bp )
  print( cp )
}
# Any notable sites? 
# Answer:
#
# Now that we are happy with no outliers, we can check for #
# correlation among predictors. Why is this important?
cor( preddf[ , prednames] )
# Are there any predictors we need to worry about?
# What correlations would be worrisome?
### end predictor check ----------------

# Now that we are satisfied with our predictor data we can #
# append it to our new closeddf:
closeddf <- #select the columns we want to keep in preddf 
  preddf %>% dplyr::select( o.sites, year, all_of(prednames) ) %>%
  right_join( closeddf, by = c("o.sites", "year") )
# why did we use right_join instead of left_join?
# check output #note that I always check dimensions when joining
# dataframes. Sometimes we add or substract rows unintentionally 
# so it is good to check that your code had the desired result
head( closeddf); dim( closeddf )

########################################################################
# We repeat the process for our robust design. we want to join #
# the two dataframes keeping relevant columns 
head( obs_df )
opendf <-  obs_df %>%
  #select desired columns to keep:
  dplyr::select( o.sites, year, pres.j1, pres.j2, pres.j3,
                 observer.j1, observer.j2, observer.j3 ) 
#check
head( opendf ); dim( opendf )
# We append predictors:
opendf <- preddf %>% 
  dplyr::select( o.sites, year, all_of(prednames) ) %>%
  left_join( opendf, by = c("o.sites", "year") )
#check
head( opendf ); dim( opendf )
# We also check for missing values in the response #
# as those are often not allowed in frequentist analyses
colSums( is.na( opendf[, c("pres.j1", "pres.j2","pres.j3")]) )
#none are present

################################################################
##########    save relevant data and workspaces     ###########
#save closed dataframe in our data folder:
write.csv( closeddf, paste( getwd(),"/Data/closedf.csv", sep = "" ),  
           row.names = FALSE )

#save open dataframe in our data folder:
write.csv( opendf, paste( getwd(),"/Data/opendf.csv", sep = "" ),  
           row.names = FALSE )

# if you want to save any of the plots we produced #
# for your presentation or ms, you do it here too:
# Save the most recently viewed plot with ggsave() to define file type, resolution, 
# and plot dimensions:
ggsave("Data/AprMayTXYear.png", dpi=500, 
       height = 10, width = 15, units= "cm" )
# or if you saved it as an object:
#start by calling the file where you will save it
tiff( 'Data/FebTXYear.tiff',
       height = 10, width = 15, units = 'cm', 
      compression = "lzw", res = 400 )
#call the plot
bp
#turn off
dev.off()

# if you want to save your workspace, because you are still #
# working through it use the following command:
#save.image( "DataPrepWorkspace.RData" )
########## End of saving section ##################################
################## Save your data and workspace ###################

############# END OF SCRIPT ########################################