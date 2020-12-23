#######################################################################
#######################################################################
####   This script was created by Dr. Jen Cruz as part of      ########
####          the Applied Population Ecology Class              #######
####                                                           ########
#### In this script we simulate the data that we will use      ########
#### for all the sample scripts provided in this course.       #######
#### We simulate three datasets: occupancy records for multiple ######
#### sites over multiple years. This will involve repeat surveys #####
###  within a primary season when the populations are assumed to #####
###  be closed, with the population assumed to be open between  #####
#### primary seasons. This is the classic ROBUST DESIGN.         #####
### We will use a robust design to also simulate sampling over a #####
### subset of these sites using repeated counts. Lastly, another ####
### smaller subset of sites will be simulated to be sampled using ####
### capture-recapture methods, again following a robust design.   #####
####                                                             ##### 
### We link these demographic parameters (occupancy, abundance,   ###
### site persistence and extinction, survival) to predictors    ####
### We also link detection to predictors.                         ####
#####################################################################
##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

# Install new packages from "CRAN" repository. # 
# If you previously installed the packages, load the  ones you need # 
# using the library() function.You only need to install packages once, #
# however, you must reload the packages every time you start up RStudio. # 

# install packages 
install.packages( "tidyverse" ) #actually a collection of packages 
# including dplyr, lubridate, tidyr, ggplot2 and more.

#install.packages( "unmarked" )
install.packages( "sf" )
install.packages( "rgdal")
install.packages( "raster" )

# load packages:
library( tidyverse ) 
#library( unmarked )
library( sf )
library( sp )
library( rgdal )
library( raster )
## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- "C:/Users/jencruz/Google Drive/QCLabShared/Data/"

# load NCA polygon
NCA <-  st_read( paste( datadir, "NCA/GIS_NCA_IDARNGpgsSampling/BOPNCA_Boundary.shp", 
        sep = "" ), quiet = TRUE )
#quick plot
ggplot( NCA ) + geom_sf( aes( geometry = geometry ) )

# load climate data from Twin Falls station
climraw <- read.csv( file = paste( datadir, "Climate/TwinFallsClimateData.csv", sep = ""),
                     header = TRUE )
#check
head( climraw )

#import sagebrush layer 
# start by setting pathway
sagepath <- paste( datadir, "Habitat/Sage_2007_2018/", sep = "" )
# extract file names
sagefiles <- dir( sagepath, pattern = "rec_v1.img", full.names = TRUE  )
sagefiles
#get names only
sagenames <- dir( sagepath, pattern = "rec_v1.img", full.names = FALSE  )
sagenames

#import annual grasses layer:
# start by setting pathway
grasspath <- paste( datadir, "Habitat/Grass_2007_2018/", sep = "" )
# extract file names
grassfiles <- dir( grasspath, pattern = "rec_v1.img", full.names = TRUE  )
grassfiles
#get names only
grassnames <- dir( grasspath, pattern = "rec_v1.img", full.names = FALSE  )
grassnames
########## end of data load ################################
#######################################################################
#### Simulating occupancy data #################
# Our study species is the Piute ground squirrel, Spermophilus    # 
# mollis. Ground squirrels are widely distributed in sagebrush steppe #
# habitats of the Great Basin and Columbia Plateau. #
# Their abundance is influenced by drought, low temperatures when #
# they emerge from hibernation in Feb and high temperatures in April-May #

# Think of this simulation as designing your field study. We first #
# need to define how many sites to sample, how many repeat surveys #
# each primary season (i.e.,season when the population is assumed #
# to be closed) and how many primary seasons (1 or more)?         #
# Let's start defining these:
# Number of years we will run the study (primary seasons):
# study will run from 2007 to 2018. Define year range:
yrrange <- 2007:2018
# create column names for each year
yrnames <- paste( "yr", yrrange, sep = "")
#number of primary seasons:
T <- length( yrrange )
# Number of repeat surveys each season:
J <- 3
#Note that at least 3 surveys are recommended when detection is >0.5
# Mackenzie and Royle (2005) J. App. Ecol. 42: 1105-1114.

# Number of sites we visit to search for presence:
Io <- 100
# Number of sites we visit to count individuals observed, without marking them:
Ic <- 50
# Number of sites we visit to trap and mark individuals:
Im <- 20

# Start by randomly selecting sites to sample:
# choose sites inside NCA, note some may be too close to be independent so #
# we start with a bigger selection than what we need
sites <- st_sample( NCA, size = Io*2, exact = TRUE )
#calculate distances between sites
site.dist <- st_distance( sites )
# replace 0 with NA to calculate if any are closer than 500m
diag( site.dist) <- NA
# remove those that are closer than 500 m:
keep <- which( apply( X = site.dist, MARGIN = 1, FUN = min, na.rm = TRUE ) > 500 )
#from this reduced sample then randomly select those for occupancy sampling 
keep <- sample( keep, size = Io )
#create a spatial point object so that we can plot them:
sitesIo <- sites[ keep, ]
#now subsample sites for counts
keepc <- sample( keep, size = Ic )
sitesIc <- sites[ keepc, ]
# now subsample sites for trapping
keepm <- sample( keepc, size = Im )
sitesIm <- sites[ keepm, ]
#check:
check <- head( site.dist )
diag(check) <- NA
apply( X = check, MARGIN = 1, FUN = min, na.rm = TRUE )

#plot site locations:
ggplot( NCA, aes( geometry = geometry ) ) + 
  theme_bw( base_size = 15 ) +
  geom_sf( size = 1.5) +  geom_sf( data = sites, size = 1 ) +
  geom_sf( data = sitesIo, size = 1.5, color = "blue" ) +
  geom_sf( data = sitesIc, size = 2, color = "green" ) +
  geom_sf( data = sitesIm, size = 3, color = "orange" )

# Create an sf dataframe that will store our occupancy data by combining #
#our simple features (spatial points) with original site IDs:
occdf <-  st_sf( geometry = st_sfc( sitesIo ), orgID = keep )
#view
head( occdf)
#add new siteID for occupancy sites:
occdf[ ,'o.sites' ] <- 1:length( keep )
#view
head( occdf )
# add columns as to whether they were sampled with point counts
occdf[ ,'counted' ] <- 'no'
#yes if they were sampled with point counts
occdf[ occdf$orgID %in% keepc,'counted' ] <- 'yes'
#add columns for trapped sites
occdf[ ,'marked' ] <- 'no'
#yes if they were trapped
occdf[ occdf$orgID %in% keepm,'marked' ] <- 'yes'
#check
head( occdf )

#work out new site ID for count sites:
c.sites <- occdf$o.sites[ occdf$orgID %in% keepc ]
#work out new site ID for mark sites:
m.sites <- occdf$o.sites[ occdf$orgID %in% keepm ]

#create sf dataframe to store our point count data:
countdf <- st_sf( c.sites = c.sites, geometry = st_sfc( sitesIc ) )
#create sf dataframe to store our trapping data:
markdf <- st_sf( m.sites = m.sites, geometry = st_sfc( sitesIm ) )
#check
head( countdf ); head( markdf )

#########get climate data ####
#view raw data:
head( climraw )
str( climraw )
# We use lubridate first to get date in the correct format
climraw$prettydate <- lubridate::as_date( climraw$DATE )
# Then to extract year
climraw$year <- lubridate::year( climraw$prettydate )
# Extract month:
climraw$month <- lubridate::month( climraw$prettydate, label = TRUE, abbr = TRUE )
#view
head( climraw )
# now we summarize the data we need. Note this is assumed to be the same #
# among sites:
climdf <- climraw %>%
  #group by month
  group_by( year, month ) %>%
  # summarise minimum and maximum temperature:
  summarise( minT = min( TMIN, na.rm = TRUE ),
             maxT = max( TMAX, na.rm = TRUE ) )
#view
tail( climdf )
#select only the months we are interested in
# note we specify the base package here because the raster package #
# also has a subset function:
climdf <- base::subset( climdf, month %in% c("Feb", "Apr", "May" ) )
# now convert to long format
climdf.long <- climdf %>% gather( "minT", "maxT", key = predictor, value = value )
tail( climdf.long )
#calculate max temperature over Apr-May for each year:
climdf <- climdf.long %>% group_by( year ) %>% 
  subset( month != "Feb" & predictor == 'maxT' ) %>%
  summarise( AprMay.maxT = max( value ) )
#view 
head( climdf )
#select min temperature for Feb
minT <- climdf.long %>% group_by( year ) %>% 
  subset( month == "Feb" & predictor == 'minT' ) %>%
  summarise( Feb.minT = value )
#view
head( minT )
#join with other predictor
climdf <- left_join( climdf, minT, by = "year" )
#view 
head( climdf )
#plot 
ggplot( climdf ) + theme_bw( base_size = 15 ) +
#  geom_histogram( aes(Feb.minT ) )
  geom_histogram( aes(AprMay.maxT ) )
#### end clim data manipulation ####

########### get habitat data #########
# extract sagebrush and annual grasses 
#start by creating dataframe to store sagebrush values
sagedf <- occdf %>% st_drop_geometry() %>%
        dplyr::select( o.sites, counted, marked  )
#add site id
sagedf[ , yrnames] <- NA
#view
head( sagedf )
sagenames
# Create a function to extract proportion of habitat from raster files that #
# are combined into a stack across years:
hab_extract <- function( df, files, site.locs, buf ){
  #required inputs:
  # df is dataframe where we will populate habitat values
  # files are the file names of the rasters for the habitat type for each year
  # site.locs are the site locations as SpatialPoints
  # buf is the buffer radius in meters over which we extract habitat
  
  #import raster files
  temp.rast <- lapply( files, raster ) #lapply( files, raster )
  #turn into stack
  temp.stack <- stack( temp.rast )
  # #convert site locations to spatial points:
  # convert site location coords to match habitat raster
  locs <- spTransform( site.locs, proj4string( temp.stack ) )
  #extract nlcd landcover for cell associated with each site, each year
  cover <- raster::extract( temp.stack, locs, buffer = buf )
  
  #get mean proportion of landcover for each site (each list), for all years (columns)
  for( i in 1:dim(locs@coords)[1] ){ #loop over each site
  df[i,yrnames] <- apply( cover[[i]], MARGIN = 2, mean ) 
  }
  # Delete any intermediate files created by raster.
  raster::removeTmpFiles()
  
  # output dataframe:
  return( df )
}
#convert site location to spatial points
site.locs <- as_Spatial( sitesIo ) #as_Spatial( site.locs )
# populate  sagebrush dataframe
sagedf <- hab_extract( df = sagedf, files = sagefiles, 
                       site.locs = site.locs, buf = 200 )


#### end of habitat prep ###########

#############end of section creating data #########################
################## Save your data and workspace ###################

########## End of saving section ##################################
###################   END OF SCRIPT  ################################