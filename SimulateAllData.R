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

#install.packages( "unmarked" ) #occupancy and abundance estimates
install.packages( "sf" ) #easy spatial dataframes
install.packages( "rgdal") #import spatial data
install.packages( "raster" ) #import, extract, manipulate rasters
install.packages( "rasterVis" ) # plotting rasters
install.packages( "RColorBrewer" ) #plotting colors
install.packages( "wiqid" ) #quick estimates for wildlife populations 
# load packages:
library( tidyverse ) 
#library( unmarked )
library( sf )
library( sp ) #commonly used spatial functions
library( rgdal )
library( raster )
library( rasterVis )
library( RColorBrewer )
library( wiqid )

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
#### Simulating data #################
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
#create a dataframe of years sampled
yrdf <- data.frame( year = yrrange, yearname = yrnames )

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
climraw[ which( climraw$year == 2014 & climraw$month == 'Jun'),]
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
#select min temperature for Feb
minT <- climdf.long %>% group_by( year ) %>% 
  subset( month == "Feb" & predictor == 'minT' ) %>%
  summarise( Feb.minT = value )
#view
head( minT, 8 )
#calculate max temperature over Apr-May for each year:
maxT <- climdf.long %>% group_by( year ) %>% 
  subset( month != "Feb" & predictor == 'maxT' ) %>%
  summarise( AprMay.maxT = max( value ) )
#view 
head( maxT, 8 )

#join with other predictor
climdf <- left_join( minT, maxT, by = "year" )
#view 
head( climdf,8 )
#data are missing for Apr-May maxT we will use mean for 2013, 2015 instead
climdf$AprMay.maxT[ climdf$year == 2014 ] <- mean( climdf$AprMay.maxT[ climdf$year == 2013 ],
                            climdf$AprMay.maxT[ climdf$year == 2015 ] )
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

### combine and standardize predictor data ####
head( sagedf )
#start with habitat: converting to long format:
preddf <- sagedf %>% 
  gather( key = yearname, sagebrush, yrnames ) %>%
  left_join( yrdf, by = "yearname" ) 
#check
head( preddf ); dim( preddf )
#add cheatgrass
head( grassdf )
preddf <- grassdf %>% select( -counted, -marked ) %>%
  gather( key = yearname, cheatgrass, yrnames ) %>%
  left_join( preddf, by = c("yearname", "o.sites" ) )
#check
head( preddf ); dim( preddf )
#add climate
head( climdf )
preddf <- left_join( preddf, climdf, by = "year" ) 
#check
head( preddf ); dim( preddf )
#define predictor columns for demographic models:
prednames <- c( "cheatgrass", "sagebrush", "Feb.minT","AprMay.maxT" )
#reorder columns
preddf <- preddf %>% dplyr::select( o.sites, counted, marked, yearname, year,
                      cheatgrass, sagebrush,  Feb.minT, AprMay.maxT )

#check correlation among predictors
round( cor( preddf[ ,prednames ] ), 1)

#standardize predictors so that coefficient estimates are comparable:
#climate data only for those years we sample squirrels:
clim.std <- climdf[ 1:T, ]
#standardize each predictor column:
clim.std[,"Feb.minT"] <- wiqid::standardize( pull(clim.std, Feb.minT) )
clim.std[,"AprMay.maxT"] <- wiqid::standardize( pull(clim.std, AprMay.maxT) )
tail( clim.std)

#habitat dataframes
grass.std <- grassdf
grass.std[ ,yrnames] <- wiqid::standardize( as.matrix( grass.std[ ,yrnames] ) )
tail( grass.std )

sage.std <- sagedf
#standardize across the matrix
sage.std[ ,yrnames] <- wiqid::standardize( as.matrix( sage.std[ ,yrnames] ) )
head( sage.std )

#all predictors 
pred.std <- preddf[,c( "o.sites", "counted", "marked", "yearname", 
                      "cheatgrass", "sagebrush" ) ]
#add standardized climate data
head( pred.std )
pred.std <- left_join( pred.std, clim.std, by = "year" ) 
#standardize habitat variables
pred.std[ ,"cheatgrass" ] <- scale(pred.std[ ,"cheatgrass" ] )
pred.std[ ,"sagebrush" ] <- scale(pred.std[ ,"sagebrush" ] )
#view
head( pred.std )

#### end combine #####
######################################################################
#### simulating demographic parameters #######
head( occdf )
#create dataframes to store true occupancy and abundance values:
Odf <- Ndf <- matrix( data = 0, nrow = Io, ncol = T )
# create dataframes to store number of births and immigrants, starting with #
# the second season
B <- M <- matrix( data = 0, nrow = Io, ncol = T-1 )

# Some of our sites are going to be empty but can be colonized in the future. #
# Cheatgrass is an invasive species that increases the likelihood of # 
# more frequent fires in the system. Here we assume that more cheatgrass # 
# will increase the probability of the population becoming extinct #
# For our first season, we assumed that sites with cheatgrass >15 % start #
# out empty of ground squirrels
Odf[,1] <- ifelse( grassdf[ ,yrnames[1]] > 15, 0, 1)
# let's randomly draw abundances for the 1st year for those occupied sites:
Ndf[,1] <- round( runif( n = Io, 1, 500 )) * Odf[,1]
# We then define the dynamic occupancy process for the following seasons:
# If a site is occupied, O, it has a high probability of remaining #
# occupied, phi, except if cheatgrass increases too much, or if the plague  # 
# comes through the site.

# Define the relationship between the probability of remaining occupied, phi, 
# and cheatgrass
#intercept: 
int.phi <- 1
# coefficient for cheatgrass
beta.phi <- -4
# relationship with phi with a logit link
logit.phi <- int.phi + ( beta.phi * grass.std[ ,yrnames ] )
# let's plot the relationship for one year to see what it looks like:
phidf <- data.frame( ID =  1:Io, 
                     # use non-standardized values
                     x = grassdf[ ,yrnames[t] ], 
                     # we calculate inverse logit using plogis:
                     y = plogis( logit.phi[,yrnames[t]] )  )
#plot
ggplot( phidf, aes( x = x, y = y ) ) +
  labs( x = "Cheatgrass cover (%)", y = "Probability of remaining occupied")+
  theme_bw(base_size = 15) + geom_point( size = 2 )
# Let's assume the probability of the plague coming through is constant through #
# time and small:
plague <- 0.02

# If a site is unoccupied (1-O), its probability of becoming recolonized, gamma,
#is related to an increase in sagebrush #
#Let's define that relationship:
#Intercept
int.gamma <- -3
#coefficient for sagebrush
beta.gamma <- 2
#relationship with colonization probability with a logit link:
logit.gamma <- int.gamma + ( beta.gamma * sage.std[ ,yrnames ] )
# let's plot the relationship to see what it looks like:
gammadf <- data.frame( ID =  1:Io, x = sagedf[, yrnames[t]], 
                     y = plogis( logit.gamma[, yrnames[t]] ) )
#plot
ggplot( gammadf, aes( x = x, y = y) ) +
  labs( x = "Sagebrush cover (%)", y = "Probability of colonization") +
  theme_bw(base_size = 15) + geom_point( size = 2 )

# Now let's simulate population growth for the following years using a #
# Gompertz model adapted to discrete time steps. #
# See: Cruz et al. 2013 PLOS ONE 8(9):e73544 for example.

# Female Piute ground squirrels give birth to an average of 5-10 young #
# Reproduction is affected by food availability early in the #
# season when they come out of hibernation, with colder Feb temperatures #
# signifying less food, lower reproduction and also lower survival of #
# adults. 
# Survival is also affected by really hot temperatures, with individuals #
# unable to forage when temperatures are too hot. So we expect a #
# negative relationship between survival and max T in Apr-May #
#let's define these relationships
# Lastly, survival is expected to be higher in sites with more sagebrush #
int.psi <- log(1.1) #this intercept is the log( mean growth rate )
#coefficients for sagebrush, feb.minT, aprmay.maxT:
beta.psi <- c(  0.1, 0.4, -0.3 )  
#create coefficient vector
coefs.psi <-  as.vector( c( int.psi, beta.psi ) )
#define predictor matrix
psi.preds <- as.matrix( cbind( rep(1, dim(pred.std)[1] ),
          pred.std[ ,c("sagebrush", "Feb.minT", "AprMay.maxT") ] ) )
# matrix multiply coefficients by predictors to estimate logit.psi:
log.psi <- psi.preds %*% coefs.psi

# let's plot partial relationships to see what they look like:
for( p in 2:length( coefs.psi) ){
  psidf <- data.frame( x = preddf[ ,prednames[p] ],
                     #calculate partial prediction values
    y = exp( int.psi + ( coefs.psi[p] *     
                      as.vector( pred.std[ ,prednames[p] ] ) ) ) )
#plot partial prediction plots
a <- ggplot( psidf, aes( x = x, y = y ) ) +
  labs( x = prednames[p], y = "Population growth rate" ) +
  theme_bw(base_size = 15) + geom_point( size = 2 )
print(a )
}

# Turn growth rate into a siteXyear matrix:
log.psi.df <- cbind( preddf[ ,c("o.sites", "yearname") ], log.psi )
log.psi.df <- spread( data = log.psi.df, key = yearname, value = log.psi )
head( log.psi.df )  
#view historgram of growth rates:
hist(exp(log.psi))

# Estimate occupancy and abundance for the following time periods in a loop:
for( t in 1:(T-1) ){
  # Now we can derive occupancy as the product of the colonization and extinction #
  # processes derived above:
  #We calculate probability of remaining occupied, phi, by subtracting the 
  # plague probability:
  Phi <- plogis( logit.phi[,t] ) - plague
  Phi <- ifelse( Phi < 0, 0, Phi )
  # We estimate colonization probability gamma:
  Gamma <- plogis( logit.gamma[,t] )
  # we derive occupancy 
  Odf[ ,t+1] <- rbinom( n = Io, size = 1, prob = ( ( Odf[,t] * Phi ) + 
                ( (1 - Odf[,t]) * Gamma ) ) )

  # Abundance, N[t+1] is determined by a Gompertz process driven by site #
  # occupancy, O[t+1], abundance in the previous season, N[t], population #
  # growth rate, Psi[t] and density-dependence, with Poisson error. #
  # In addition, potential migrants may be added to a site depending on 
  # a binomial process driven by a random draw of dispersers and the probability #
  # that the site was colonized that year, Gamma[t]. 

  #Estimate the population growth rate for that site, that year, Psi, by adding
  # a density-dependent term when the site was occupied in the previous season:
  Psi <- exp( log.psi.df[ ,yrnames[t]] + ( -0.05 * log( Ndf[,t] + 1 ) ) )

  # Determine the number of survivors, based on current site occupancy, Odf[t+1], #
  # previous abundance, N[t], and population growth rate, Psi[t]. The later takes #
  # into account covariates and reflects births and deaths in the population:
  S <- rpois( n = Io, lambda = Ndf[,t] * Psi * Odf[,t+1] )

  # Define potential migrants, M, as the product of a binomial process drawing #
  # from a random draw of potential dispersers ranging from 1 to 20, and #
  # the probability of colonization for that year for each site: #
  M <- rbinom( Io, size = round(runif( Io, min = 1, max = 20 )), prob = Gamma )

  # Calculate abundance as the sum of survivors, S, and migrants, M:
  Ndf[ ,t+1 ] <- S + M 

} #end of loop

#check output
Ndf
head( Odf )
rowSums( Odf )
#######
#### turn matrices to sf dataframes to save them and to long formats: ######
# True occupancy dataframe
Odf.save <- as.data.frame( Odf )
#add column names
colnames( Odf.save ) <- yrnames
#view
head( Odf.save )
#add site attributes
Odf.save <- cbind( occdf, Odf.save )

# True abundance dataframe:
Ndf.save <- as.data.frame( Ndf )
#add column names
colnames( Ndf.save ) <- yrnames
#view
head( Ndf.save )
#add site attributes
Ndf.save <- cbind( occdf, Ndf.save )

# reformat true occupancy and abundance dataframes to a long format:
head( Odf.save )
# extract columns from Odf and drop geometry:
Odf.long <- st_drop_geometry( Odf.save[ ,c('o.sites',
                                           'counted', 'marked', yrnames ) ] )
#convert to long formats
Odf.long <- Odf.long %>% 
  gather( key = yearname, Occ, yrnames )
#check that it worked
#view
head( Odf.long ); dim( Odf.long )

# repeat process for true abundance dataframe
Ndf.long <- st_drop_geometry( Ndf.save[ ,c('o.sites','counted', 'marked', yrnames ) ] )
#convert to long formats
Ndf.long <- Ndf.long %>% 
  gather( key = yearname, N, yrnames )
#check that it worked
#view
head( Ndf.long ); dim( Ndf.long )
#### end matrix to df conversion ##########

########################################################################
##### add imperfect detection to our sampling ######
# Create an observation dataframe to store results:
#survey names
jnames <- c( "j1", "j2", "j3" )
# create column names where we will store observations
onames <- c( paste( "pres", jnames, sep = "." ),
             paste( "count", jnames, sep = "." ) )
#create dataframe with columns that will contain observations for the three #
# sampling methods
obsdf <- data.frame( matrix( 0, nrow = dim(preddf)[1], ncol = length( onames ) ,
                             dimnames = list(NULL, onames ) ) )
# add ID information 
obsdf <- cbind( preddf[ ,c('o.sites', 'counted', 'marked', 'yearname', 'year')],
                obsdf )
### simulate detection predictors at the survey, J, level or animal, A, level ####
# observer effects:
#let's create ID for 4 observers
observers <- paste( "tech", 1:4, sep = "." )
#create column ids for observer predictors
obscols <- paste( 'observer', jnames, sep = "." )
# create column id for time of day:
timecols <- paste( "time", jnames, sep = "." )
#create range of times in minutes ranging from 0 at 6am to 360 at noon:
survtimes <- seq( 0, 360, 5 )

#loop to create random choice of observer for each survey, site, season:
for( j in 1:J ){
obsdf[ ,obscols[j] ] <- sample( x = rep(observers, dim( obsdf)[1]) ,
                               size = dim( obsdf)[1], replace = FALSE )
}
#loop over surveys to create random time that they were conducted:
for( j in 1:J ){
obsdf[ ,timecols[j] ] <- sample( x = survtimes, dim( obsdf)[1], replace = TRUE )
}

#check
head( obsdf )

### end detection predictors ######
#### Detection for occupancy searches ######
# detection for occupancy searches is related to % sagebrush, with  #
# increasing sagebrush lowering visibility of ground squirrel sign #
# and to observer effects with 4 observers randomly surveying sites #
# on different days #

# Let's define the relationship between sagebrush and detection:
#intercept as the logit of mean detection:
int.p.occ <- qlogis( 0.7 )
#coefficient for sagebrush[i,t]
beta.p.occ <- -0.3
#combine
logit.p <- int.p.occ + ( beta.p.occ * pred.std[, 'sagebrush'] ) 
# let's plot the relationship to see what it looks like:
#combine predictor dataframe (not standardized) with response:
cbind( preddf, logit.p ) %>% 
  #add y as expit(logit.p )
  mutate( y = plogis( logit.p ) ) %>%
  #plot result
ggplot( ., aes( x = sagebrush, y = y) ) +
  labs( x = "Sagebrush cover (%)", y = "Probability of detection") +
  theme_bw(base_size = 15) + geom_point( size = 2 )

hist( preddf$sagebrush)
#now let's add the observer effect:
obsefs <- data.frame( id = observers, effect = rnorm(n = 4, mean = 0, sd = 1 ) )
#calculate the mean probability of detection for each observer
plogis( int.p.occ + obsefs$effect )
obsefs
# we can now populate the observations:
for( j in 1:J ){
  for( i in 1:dim(obsdf)[1]){
    # calculate detection for that survey by adding the respective observer effect #
    # to the logit.p relationship with sagebrush
    p.j[i] <- plogis( logit.p[i] + 
        obsefs[ which( obsefs[,"id"] == obsdf[ i, obscols[j] ] ), "effect" ] )

    # now populate observed presence as a Binomial process based on the #
    # detection probability, p[i,t], and true occupancy, Occ[i,t] #
    obsdf[ i,onames[j] ] <- rbinom( n = 1, size = 1, 
                                 prob = p.j[i] * Odf.long[ i,"Occ" ] )
  }
  # print(p.j[1:100])
  # ap <- hist( p.j )
  # bp <- hist( p.j * Odf.long[,"Occ"] )
  # cp <- table( obsdf[,onames[j]] )
  # print(cp)
}
#view
tail( obsdf,10 )
colSums( obsdf[,onames[1:3]] )
#### end detection for occupancy ####
#### Detection for Point Counts ######
# detection for point counts is related to time of day #
#start by creating standardized matrix of our survey times:
times.std <- wiqid::standardize( as.matrix( obsdf[ ,timecols] ) )
#view
head( times.std )
#define our relationship with time of day as a quadratic:
#intercept is logit of mean detection
int.p.count <- qlogis( 0.7 )
#coefficients for quadratric relationship with time of day
beta.p.count <- c( -0.2, -0.4 )
logit.p.count <- int.p.count + ( beta.p.count[1] * times.std[, timecols] ) +
                 ( beta.p.count[2] * times.std[ ,timecols ]^2 )
#view
head(logit.p.count)
#change colnames
colnames( logit.p.count ) <- paste( "logit.p", jnames, sep = ".")
# let's plot the relationship to see what it looks like:
for( j in 1:J){
  #combine predictor dataframe (not standardized) with response:
jp <-  data.frame( x = obsdf[ ,timecols[j] ], 
              y = plogis( logit.p.count[,j] ) ) %>% 
  #plot result
  ggplot( ., aes( x = x, y = y) ) +
  labs( x = "Survey time", y = "Probability of detection", main = jnames[j]) +
  theme_bw(base_size = 15) + geom_point( size = 2 )
print( jp )
}

dim(obsdf );dim(Ndf.long); dim( logit.p.count)
# we can now populate count observations:
for( j in 1:J ){
  for( i in 1:dim(obsdf)[1]){
    # observed counts are a Poisson process dependent on detection probability #
    # and true abundance:
    obsdf[ i,onames[3+j] ] <- rpois( n = 1, lambda = plogis( logit.p.count[i,j] ) * 
                                  Ndf.long[ i,"N" ] )
  }}
head( obsdf,30 )
######
#### Detection for Trapping ######
# detection for trapping is related to trap happiness and sex #
# let's assume and even sex ratio and randomly assign individuals as male #
# or female
#define mean probability of 1st capture for female:
int.p.mark <- qlogis( 0.6 )
#define probability of recapture and of capture if male:
beta.p.mark <- c( 1, 0.5 )
p.mark <- data.frame( sex = rep( c( "female", "male"),2 ), trap.resp = c("no", "no","yes","yes") )
#populate with detection:
#detection for 1st capture, female
p.mark[1,"p"] <- plogis( int.p.mark )
#relationship for 1st capture male
p.mark[2,"p"] <- plogis( int.p.mark + beta.p.mark[2] )
#relationship for recapture female
p.mark[3,"p"] <- plogis( int.p.mark + beta.p.mark[1] )
#detection for recapture, male
p.mark[4,"p"] <- plogis( int.p.mark + sum(beta.p.mark ) )
p.mark

#let's create a dataframe of true abundance by selecting sites that were trapped:
Ndf.mark <- dplyr::filter( Ndf.long, marked == "yes" ) %>%
        dplyr::select( -counted )
head( Ndf.mark )

#create dataframe to store data for individuals available for capture
Indf <- data.frame( o.sites = rep(Ndf.mark[1,'o.sites'] , Ndf.mark$N[1] ), 
          year = rep( yrdf$year[ which( yrdf[ ,'yearname'] == Ndf.mark[1,'yearname'] ) ], 
                                Ndf.mark$N[1] ) )
#assign gender
Indf[ ,"sex"] <- sample( c("female", "male"), size = Ndf.mark[ 1,"N" ], 
                                        replace = TRUE, prob = c(0.5, 0.5 ) )
#view
head( Indf );dim(Indf)

#create a row for every individual available for capture:
for( i in 2:dim(Ndf.mark)[1] ){
  if( Ndf.mark$N[i] > 0 ){
    #assign gender
    sex <- sample( c("female", "male"), size = Ndf.mark$N[i], 
                           replace = TRUE, prob = c(0.5, 0.5 ) )
    #extract site id:
    o.sites <- Ndf.mark$o.sites[i] 
    #extract year id
    year <- yrdf$year[ which( yrdf$yearname == Ndf.mark$yearname[i] )]
    #create temp dataframe
    df <- data.frame( o.sites, year, sex )
    #add rows to Individual data frame:
    Indf <<- rbind( Indf, df )
  } #close if function
  }# close for loop
#view
head( Indf ); dim( Indf )

#now iterate through individuals to create their capture history
trapnames <- paste( "trap", jnames, sep = "." )
for( j in 1:3 ){
  Indf[ ,trapnames[j] ] <- 0
}
head( Indf )
#populate trap history for all individuals that were available for capture:
for( i in 1:dim(Indf)[1]){
  for( j in 1:J ){
    #select trap response based on whether individuals was previously captured that 
    #season. Note that they are all new captures to the season  on the 1st survey#
    # so there is no behavioral response for j1:
    resp <- ifelse( sum( Indf[ i,trapnames ] ) > 0, "yes", "no" )
    #choose detection probability based on sex and previous capture of individual
    p.i <- p.mark$p[ which( (p.mark$sex == Indf$sex[i]) & (p.mark$trap.resp == resp) ) ]
    #use binomial draw to work out if it was trapped
    Indf[i, trapnames[j]] <- rbinom(n = 1, size = 1, prob = p.i   )
}}
#now add capture history:
Indf$ch <- apply( Indf[,trapnames ], 1, paste, collapse="" )
#check
head( Indf ); dim( Indf )
#remove those individuals that were never captured:
Indf <- filter( Indf, ch != "000" )
#add individual ids:
Indf$indid <- 1:dim(Indf)[1]
#check
head( Indf ); dim( Indf )

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
#save shapefile of site locations
sf::st_write( occdf, paste( getwd(), "/Data/sites.shp",  sep = "" ),  
              driver = "ESRI Shapefile" )

#save predictor dataframe in longformat
write.csv( preddf, paste( getwd(),"/Data/predictors.csv", sep = "" ),  
                          row.names = FALSE )

#save standardized predictors
write.csv( pred.std, paste( getwd(),"/Data/predictors_std.csv", sep = "" ),  
           row.names = FALSE )

#save true occupancy dataframe as shapefile so that we can keep spatial information:
sf::st_write( Odf.save, paste( getwd(), "/Data/Odf.shp", sep = "" ),
              driver = "ESRI Shapefile" )

#save true abundance dataframe as shapefile so that we can keep spatial information:
sf::st_write( Ndf.save, paste( getwd(), "/Data/Ndf.shp", sep = "" ),
              driver = "ESRI Shapefile" )

#save presence and count observations:
write.csv( obsdf, paste( getwd(),"/Data/obsdf.csv", sep = "" ),  
           row.names = FALSE )

#save individual captures
write.csv( Indf, paste( getwd(),"/Data/Indf.csv", sep = "" ),  
           row.names = FALSE )

#save workspace 
save.image( "SimDataWorkspace.RData" )
########## End of saving section ##################################
###################   END OF SCRIPT  ################################