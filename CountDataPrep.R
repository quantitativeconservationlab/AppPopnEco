#######################################################################
#######################################################################
####   This script was created by Dr. Jen Cruz as part of      ########
####          the Applied Population Ecology Class              #######
####                                                           ########
# Our study species is the Piute ground squirrel, Spermophilus        # 
# mollis at the Morley Nelson Birds of Prey Nature Conservation Area. #
# Their abundance is influenced by drought, low temperatures when     #
# they emerge from hibernation in Feb and high temperatures in        #
# in April-May                                                        #
# 50 sites were randomly selected for repeated count surveys over     #
# three days. Surveys were repeated over multiple years ()            #
# Surveys involved point counts for 2min where all individuals detected #
# over a 200 m radius were recorded.                                  # 
#######################################################################
##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

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

# We start with our robust design data ##########
# Remember that not all sites were surveyed with point counts but only # 
# those were 'counted' column = yes. We start by filtering the dataframe #
# to sites that were counted. We also remove presence/absence records. #
# The same observers that conducted the presence surveys, conducted #
# the counts so those columns remain. Note that time of day can now be #
# included as they reflect the time that the point counts were conducted #

# Also note that our siteXyear predictor data is the same as what we #
# used in the occupancy analyses. We therefore do not explore those for #
# correlation, outliers, spread, since we have already done it. # 
# Refer to OccDataPrep.R for details of how to do that if you are #
# working with your own data. #

head( obs_df )
opendf <-  obs_df %>% 
  # keep only rows where sites had point counts:
  dplyr::filter( counted == "yes" ) %>%
  #select columns to keep:
  dplyr::select( o.sites, counted, year, count.j1, count.j2, count.j3,
                 observer.j1, observer.j2, observer.j3,
                 time.j1, time.j2, time.j3 ) 
#check
head( opendf ); dim( opendf )

# We append predictors:
opendf <- preddf %>% 
  dplyr::select( o.sites, year, cheatgrass, sagebrush, Feb.minT, AprMay.maxT ) %>%
  right_join( opendf, by = c("o.sites", "year") )
#check
head( opendf ); dim( opendf )

# Why do we use right_join instead of left_join or other?
# Answer: 
#

# Let's define some parameters for our robust design:
# How many sites were sampled?
I <- length( unique( opendf$o.sites ) )
# What years were sampled (primary seasons)?
yrrange <- sort( unique( opendf$year) )
#How many years (i.e., primary seasons )?
T <- length( yrrange )
# How many replicate surveys
J <- 3
#view
I; yrrange; T; J

###########

# We now need to decide which year to use for our closed analysis ###
# let's work out what our counts look like each year:
# we use our opendf and group by year by site:
opendf %>% group_by( year, o.sites ) %>%
  # we calculate the maximum count observed during the primary season
  # this is our 'naive' estimate of abundance for that site that year:
  # why?
  mutate( maxc = max( count.j1, count.j2, count.j3 ) ) %>%
  # we select columns we need for plotting
  select( o.sites, year, maxc ) %>% ungroup() %>%
  # we group data by year
  ggplot(., aes(group = year )) +
  # choose a preset theme
  theme_bw( base_size = 15 ) + 
  # plot histogram of abundances: 
  geom_histogram( aes( maxc ), binwidth = 10 ) +
  # make a different plot for each year:
  facet_wrap( ~year, ncol = 3, scales = "free" )

#how are counts changing at each site
opendf %>% group_by( year, o.sites ) %>%
  mutate( maxc = max( count.j1, count.j2, count.j3 ) ) %>%
  dplyr::select( o.sites, year, maxc ) %>% ungroup() %>%
  # we group data by year
  ggplot(., aes(color = as.factor(o.sites), x = year, y = maxc )) +
  theme_bw( base_size = 15 ) + 
  geom_line()
  
# What do we see? What is happening to site abundances, but also to #
# site occupancy? 
# Answers:
#

# We need to choose a year for our single season (closed population) analysis.
# Based on the plots I will choose 2008. I think later years may require #
# a zero-inflated Poisson model to account for all the empty sites #

# We create our closeddf:
closeddf <- opendf %>% 
  #filter only rows for desired year:
  dplyr::filter( year == 2012 ) 

################################################################
##########    save relevant data and workspaces     ###########
#save open dataframe in our data folder. Make sure not to overwrite #
# your occupancy files:
write.csv( opendf, paste( getwd(),"/Data/open_counts.csv", sep = "" ),  
           row.names = FALSE )

#save closed dataframe in our data folder:
write.csv( closeddf, paste( getwd(),"/Data/closed_counts.csv", sep = "" ),  
           row.names = FALSE )

# if you want to save your workspace, because you are still #
# working through it use the following command:
save.image( "CountPrepWorkspace.RData" )
########## End of saving section ##################################
################## Save your data and workspace ###################

############# END OF SCRIPT ########################################