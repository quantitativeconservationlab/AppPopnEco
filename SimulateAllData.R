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

install.packages( "unmarked" )
install.packages( "sf" )
install.packages( "rgdal")
install.packages( "raster" )
install.packages("remotes")
remotes::install_github("juoe/spatialtools")

# load packages:
library( tidyverse ) 
library( unmarked )
library( sf )
library( rgdal )
library( raster )

# set datadir 
datadir <- "C:/Users/jencruz/Google Drive/QCLabShared/Data/"

getwd()  
###################################################################
#### Load or create data -----------------------------------------
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
habpath <- paste( datadir, 
"NLCD_LandCoverChange_Index/NLCD_Land_Cover_Change_Index_L48_20190424.img", sep = "")
#import tiff as raster file using raster package:
habrast <- raster::raster( habpath )
#quick check by plotting with base
plot( habrast )
habrast

#### Simulating occupancy data #################
# Our study species is the Piute ground squirrel, Spermophilus    # 
# mollis. Ground squirrels are widely distributed in sagebrush-steepe #
# habitats of the Great Basin and Columbia Plateau. #
# Their abundance is influenced by drought, low temperatures when #
# they emerge from hibernation in Feb and high temperatures in April-May #

# Think of this simulation as designing your field study. We first #
# need to define how many sites to sample, how many repeat surveys #
# each primary season (i.e.,season when the population is assumed #
# to be closed) and how many primary seasons (1 or more)?         #
# Let's start defining these:
# Number of years we will run the study (primary seasons):
T <- 10
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


# populate sagebrush
occdf <- occdf %>%
      mutate( sagebrush = runif( n = Io, min = 0, max = 50 ),
              cheatgrass = runif( n = Io, min = 0, max = 50 ), 
              nativegrass = runif( n = Io, min = 0, max = (sagebrush - cheatgrass )),
              other = 1 - (sagebrush + cheatgrass + nativegrass ) )





#############end of section creating data #########################
################## Save your data and workspace ###################

########## End of saving section ##################################
###################   END OF SCRIPT  ################################