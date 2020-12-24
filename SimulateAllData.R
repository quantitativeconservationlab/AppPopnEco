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
install.packages( "sf" ) #easy spatial dataframes
install.packages( "rgdal") #import spatial data
install.packages( "raster" ) #import, extract, manipulate rasters
install.packages( "rasterVis" ) # plotting rasters
install.packages( "RColorBrewer" ) #plotting colors
# load packages:
library( tidyverse ) 
#library( unmarked )
library( sf )
library( sp ) #commonly used spatial functions
library( rgdal )
library( raster )
library( rasterVis )
library( RColorBrewer )

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
sagefiles <- dir( sagepath, pattern = "_v1.img", full.names = TRUE  )
sagefiles
#get names only
sagenames <- dir( sagepath, pattern = "_v1.img", full.names = FALSE  )
sagenames

#import annual grasses layer:
# start by setting pathway
grasspath <- paste( datadir, "Habitat/AnnualHerb_2007_2018/", sep = "" )
# extract file names
grassfiles <- dir( grasspath, pattern = "_v1.img", full.names = TRUE  )
grassfiles
#get names only
grassnames <- dir( grasspath, pattern = "_v1.img", full.names = FALSE  )
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
# define the radius (in m) of our effective sampling area
buf <- 200

# Start by randomly selecting sites to sample:
# choose sites inside NCA, note some may be too close to be independent so #
# we start with a bigger selection than what we need
sites <- st_sample( NCA, size = Io*2, exact = TRUE )
#calculate distances between sites
site.dist <- st_distance( sites )
# replace 0 with NA to calculate if any are closer than 500m
diag( site.dist) <- NA
# remove those that are closer than 500 m:
keep.org <- which( apply( X = site.dist, MARGIN = 1, 
                FUN = min, na.rm = TRUE ) > 500 )
#now need to check that these locations are not missing habitat data:
#bring in a sample raster assuming equal coverage across all
samp.rast <- raster::raster( sagefiles[1] )
# convert site location to spatial object
sites.locs <- as_Spatial( sites[ keep.org ] )
# transform coords system to match habitat raster
keep.locs <- spTransform( sites.locs, proj4string( samp.rast ) )
# why not the other way around?

#extract annual % landcover for cells within buffer surrounding each site #
# center point 
landcover <- raster::extract( samp.rast, keep.locs, buffer = buf )
#only keep sites without missing habitat data. Missing values are marked as 
# 101 or 102:
#create logical vector to subset original site list
keep.red <- keep.org == keep.org
for( l in 1:length(landcover) ){
  keep.red[l] <- ifelse( sum(landcover[[l]] >100) > 0, FALSE, TRUE )
}
#from this reduced sample then randomly select those for occupancy sampling 
keep <- sample( keep.org[ keep.red ], size = Io )
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
cols.sol = c("Sign" = "blue", "Count" = "green", "Trap" = "orange")
ggplot( NCA, aes( geometry = geometry ) ) + 
  theme_bw( base_size = 15 ) +
  geom_sf( size = 1.5) + # geom_sf( data = sites, size = 1 ) +
  geom_sf( data = sitesIo, size = 1.5, aes(color = "Sign")) + #color = "blue" ) +
  geom_sf( data = sitesIc, size = 2, aes(color = "Count")) + #color = "green" ) +
  geom_sf( data = sitesIm, size = 3, aes(color = "Trap")) + #color = "orange" ) +
  scale_color_manual(name = "Method", values = cols.sol ) +
  theme( legend.position = "top" )


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
#start by creating dataframe to store sagebrush and annual herbaceous values
sagedf <- grassdf <- occdf %>% st_drop_geometry() %>%
        dplyr::select( o.sites, counted, marked  )
#add columns for annual habitat data
sagedf[ , yrnames] <- NA
#view
head( sagedf )
#add  columns for annual habitat data
grassdf[ , yrnames] <- NA
#view
head( grassdf )

# Create a function to extract proportion of habitat from raster files that #
# are combined into a stack across years:
hab_extract <- function( df, yrnames, files, site.locs, buf ){
  #required inputs:
  # df is dataframe where we will populate habitat values
  # yrnames are column names in df where habitat values will go
  # files are the file names of the rasters for the habitat type for each year
  # site.locs are the site locations as SpatialPoints
  # buf is the buffer radius in meters over which we extract habitat
  
  #import raster files
  temp.rast <- lapply( files, raster ) #lapply( files, raster )
  #turn into stack
  temp.stack <- stack( temp.rast )

   # convert site location coords to match habitat raster
  locs <- spTransform( site.locs, proj4string( temp.stack ) )
  #extract nlcd landcover for cell associated with each site, each year
  cover <- raster::extract( temp.stack, locs, buffer = buf )
  #get mean proportion of landcover for each site (each list), for all years (columns)
  for( i in 1:dim(locs@coords)[1] ){ #loop over each site
    #masked areas have values of 101 and outside mapping areas have values of 102
    # need to replace those with NA before calculating means
    cover[[i]][ cover[[i]] %in% c(101,102) ] <- NA
    #print max value
    print( max( cover[[i]], na.rm = TRUE ))
    # calculate mean proportion for each primary season
    a <- apply( cover[[i]], MARGIN = 2, mean, na.rm = TRUE ) 
    # remove names and round values to 1 decimal 
    a <- round( as.numeric( a ), 1 )
    #print( a )
    df[ i, yrnames ] <- a
  }
  # Delete any intermediate files created by raster.
  raster::removeTmpFiles()
  
  # output dataframe:
  return( df )
}
#convert site location to spatial points
site.locs <- as_Spatial( sitesIo ) 

# populate  sagebrush dataframe. Notice that there is no NLCD data for 2012
# so we leave it empty for now and populate it later
sagedf <- hab_extract( df = sagedf, yrnames = yrnames[yrnames != "yr2012"], 
                       files = sagefiles, site.locs = site.locs, buf = buf )
#view 
head( sagedf )
max( sagedf[ ,yrnames ], na.rm = TRUE )

# populate  annual herbaceous (i.e., mainly cheatgrass) dataframe. 
# Notice that there is no NLCD data for 2012 so we leave it empty for now
grassdf <- hab_extract( df = grassdf, yrnames = yrnames[yrnames != "yr2012"], 
                       files = grassfiles, site.locs = site.locs, buf = buf )
#view
head( grassdf )
max( grassdf[ ,yrnames ], na.rm = TRUE )

# for now we assume 2012 is the mean between 2011 and 2013:
sagedf <- sagedf %>% 
          mutate( yr2012 = rowMeans( select(., yr2011, yr2013 ) ) )

grassdf <- grassdf %>% 
  mutate( yr2012 = rowMeans( select(., yr2011, yr2013 ) ) )

#### end of habitat prep ###########
#### simulating demographic parameters
head( occdf )
#create a dataframe to store true abundance values:
Ndf <- matrix( data = 0, nrow = Io, ncol = T )
# let's randomly draw abundances for the 1st year
set.seed( 2020 )
Ndf[,1] <- runif( n = Io,  )

round( runif( n = 10, 0, 500 ))
#############end of section creating data #########################
###################################################################
###### plotting #################################
### plot rasters #########
#plot sagebrush rasters 
#import raster files
temp.rast <- lapply( sagefiles, raster ) 
#turn into raster stack
temp.stack <- stack( temp.rast )
#transform locations to same coordinate sistem
locs <- spTransform( site.locs, proj4string( temp.stack ) )
#transform NCA polygon
nca <- spTransform( as_Spatial( NCA), proj4string( temp.stack ) )
#crop raster stack
ex <- extent( nca )
#crop raster stack to extent of NCA
plot.stack <- raster::crop( temp.stack, ex )
#mask raster stack
plot.stack <- raster::mask( plot.stack, nca )
# define a color ramp of yellow to green with 3 different levels
cols <- colorRampPalette( brewer.pal(3,"YlGn") )
#plot annual habitat rasters for sagebrush within the range of values observed
# in the NCA, which were < 30%
levelplot( plot.stack, at = seq(0,30,10), col.regions=cols, 
           names.attr = yrnames[yrnames != "yr2012"] ) +
  layer(  sp.points(locs,  col = "blue" ) ) +
  layer( sp.polygons( nca, col = "black" ) )

#plot cheatgrass rasters 
#import raster files
temp.rast <- lapply( grassfiles, raster ) 
#turn into a raster stack
temp.stack <- stack( temp.rast )
#crop raster stack to extent of NCA that we define above:
plot.stack <- raster::crop( temp.stack, ex )
#mask raster to only include the NCA
plot.stack <- raster::mask( plot.stack, nca )
# define color ramp to use in the plot
cols <- colorRampPalette( brewer.pal( 5,"YlOrRd") )
#plot annual habitat rasters for cheatgrass in the range of values observed
# in the NCA, which were < 50%
levelplot( plot.stack, 
           #set min and max proportion of habitat to be displayed
           at = seq(0,50,10), col.regions=cols,  
           #label panels excluding the missing data of 2012
          names.attr = yrnames[yrnames != "yr2012"] ) +
  #add site locations for occupancy study:
  layer(  sp.points(locs,  col = "blue" ) ) +
  #add NCA polygon 
  layer( sp.polygons( nca, col = "black" ) )

#####
#### plot predictors ####
yrdf <- data.frame( year = yrrange, yearname = yrnames )
head( sagedf )

#annual changes in habitat cover:
#sagedf %>% 
grassdf %>%  
  #subset( counted == 'yes' ) %>%
  #subset( marked == 'yes' ) %>%
  gather( key = yearname, cover, yrnames ) %>%
left_join( yrdf, by = "yearname" ) %>%
ggplot(., aes( x = as.numeric(year), y = cover, 
               color = as.factor(o.sites) )) + 
  #labs( x = "Year", y = "Sagebrush (%)" ) +
  labs( x = "Year", y = "Cheatgrass (%)" ) +
  theme_bw( base_size = 15 ) +
  geom_line( size = 1.3 ) +
  theme( legend.position = "none" )

#annual climate
head( climdf )
ggplot( climdf, aes( x = year, 
                     #y = Feb.minT ) ) +
                     y = AprMay.maxT ) ) +
  theme_bw( base_size = 15 ) +
  #labs( x = "Year", y = "Feb Minimum Temperature (C)" ) +
  labs( x = "Year", y = "Apr-May Maximum Temperature (C)" ) +
  geom_line( size = 1.3 )
  
###### end of plots #############

################## Save your data and workspace ###################

########## End of saving section ##################################
###################   END OF SCRIPT  ################################