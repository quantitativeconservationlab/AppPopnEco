#######################################################################
#######################################################################
####   This script was created by Dr. Jen Cruz as part of      ########
####          the Applied Population Ecology Class              #######
####                                                           ########
# Our study species is the Piute ground squirrel, Spermophilus        # 
# mollis at the Morley Nelson Birds of Prey Nature Conservation Area. #
#                                                                     #            
# We shift to real capture-recapture data from 20 sites collected     #
# during 2021 to 2023. Note that different sites were surveyed each year#
#                                                                     #
# Detection predictors: 1-effort: number of hours the traps were open #
# each day, 2-mean tempC_st up to the hour the traps were closed, 3-mean #
# wind up to the hour the traps were closed each day                   #
#
# Predictors for abundance: only 20 sitesXyear so be mindful of how            #
# many we can reasonably include. The data this time comes from RAP:
# https://rangelands.app/products/ #
# and includes measures of cover for: shrub, perennial, annual and the #
# sum of annual and perennial = herbaceous. Cover is for the year prior# # 
#
##############################################################################################################################################
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

#import individual level capture histories 
ch <- read.csv( paste0(datadir, "trapdata21-23wide.csv" ),
                header = TRUE)
#import raw ikj predictors
ikj_df <- read.csv( paste0( datadir, "traps_ikj_df_21-23.csv"),
                   header = TRUE )
#import iXk predictors
ik_df <- read.csv( paste0( datadir, "traps_ik_df_21-23.csv"),
                   header = TRUE )

########## end of data load ################################
#######################################################################
######## explore individual capture data #############
# Start viewing our observations dataframe:
head( ch ); dim( ch )

###define model parameters:
#how many sites?
(S <- length( unique(ch$SiteID) ) )
#how many unique siteXyear combo?
I <- max( ik_df$idno )
#which years (primary seasons) were sampled:
(yrrange <- sort( unique( ch$year ) ) )
#how many years were sampled:
(T <- length( yrrange ) )
#max number of surveys at a given site
J <- 3
#number of unique individuals
N <- max( ch$AnimalID, na.rm = TRUE )
#number of individuals (including augmented group)
M <- N * 3

#calculate naive abundance for each sex:
N_naive <- ch %>% group_by( id, year, Sex ) %>% 
  summarise( captures = n()  ) %>% ungroup()
#check
head( N_naive ); dim(N_naive )
#plot naive estimates:
ggplot( N_naive, aes( x = as.factor(year), y = captures,
                      color = as.factor(id)) ) +
  geom_point( size = 2 ) +
  theme_bw( base_size = 15) +
  theme( legend.position = "none" ) +
  facet_wrap( ~Sex, nrow = 2 )
# Can we see any sex differences?
# Answer:
#

# Plot naive annual estimates of abundance
N_naive %>% group_by( id ) %>% 
  summarise( total = sum(captures),
             year = first(year)) %>% 
  ggplot(., aes( x = as.factor(year), y = total,
                 color = as.factor(id) ) ) +
  geom_point( size = 2 ) +
  theme_bw( base_size = 15 )
#What does this tell us?
# Answer:
#
#
head(ch)

#create collapsed capture history column 
ch$ch  <- apply( ch[,c('j_1', 'j_2', 'j_3') ], 
                 1, paste, collapse="" )

#turn ch to factor for multi-season:
ch$ch  <- factor( ch$ch, 
      levels=c("001", "010", "011", "100", "101", "110", "111"))

# we can also summarize capture histories
table( ch$ch )
# what do you notice? 
# Answer:

#number of observed capture histories:
CH <- length(levels(ch$ch))
#create tally of individuals captured under each capture history
y_ik <- table( ch$idno, ch$ch ) 
class(y_ik)
#turn to dataframe
y_ik <- as.data.frame.matrix(y_ik)
#add idno column so that we can combine with ik_df
y_ik$idno <- as.numeric( rownames(y_ik) )
#now add sites with non of those capture histories (no captures)
# extracted from ik_df
y_ik <- ik_df %>% 
  dplyr::select( idno ) %>% 
  left_join( y_ik, by = "idno" ) %>% 
  #remove id column since we don't need it anymore
  dplyr::select( -idno )
# replace NA with 0s
y_ik[is.na(y_ik)] <- 0 
#view
y_ik
#calculate number of individuals captured per site
n_ik <- apply( y_ik, 1, sum )

##### prepare data for data augmentation analysis #####
# instead of summarizing to capture histories we keep information 
# for each captured individual: 
#extract id for sitesXyear with captures
capsites <-  unique( ch$idno)
#now add rows for sites with all 0 captures
y <- ik_df %>% filter( !(idno  %in% capsites ) ) %>% 
  dplyr::select( idno, id, SiteID, year  )
#append rows for sites wiht no captures to ch 
y <- bind_rows( ch, y )
#augment individual dataframe with M individuals
y <- rbind( y, matrix(NA, nrow = M-N, ncol = NCOL(ch), 
                      dimnames = list(NULL, colnames(ch))))

#need to replace NAs with all zero captures
y$j_1[ is.na( y$j_1 ) ] <- 0
y$j_2[ is.na( y$j_2 ) ] <- 0
y$j_3[ is.na( y$j_3 ) ] <- 0

tail(y)
#check sites got added
unique(y$idno)
#note that the trapped sites with zero captures are added at the top
#recalculate M since we added extra rows for sites with no captures
M <- dim(y)[1] #- N #M - N = number of augmented individuals

############### prepare predictor data ##############
#define predictors for detection
preds <- c("effort_mins", "tempC_st", "wind_kmph_st" )

#view histograms and correlation
for( n in preds){
  hist( ikj_df[,n], main = n )
}
#check correlation
cor(ikj_df[,preds])
# what do these two check tell us? are we ok moving forward?
# answer:
#

#Now we are ready to scale predictors
ij_sc <- ikj_df
#scale
for( i in preds ){
  ij_sc[ ,i ] <- scale( ij_sc[ ,i ] ) 
}
#check
head( ij_sc )
#turn to wide format
ij_wide <- ij_sc %>% select( -X ) %>% 
  pivot_wider( names_from = Survey,
     values_from = c( jday, effort_mins, tempC_st, wind_kmph_st),
     values_fill = 0 ) %>% 
  dplyr::arrange( SiteID, year )

# extract column ids for each detection predictor
widx <- grep( "wind", colnames( ij_wide ), value = FALSE)
tidx <- grep( "temp", colnames( ij_wide ), value = FALSE)
eidx <- grep( "effort", colnames( ij_wide ), value = FALSE)

head( ij_wide)

###for abundance predictors now plot histograms 
# and calculate correlation between predictors ###
## here:
#
# From this comment on which predictors you are choosing to use 
# moving forward? 
# Answer:
#

#now we are ready to scale predictors for abundance
ik_sc <- apply( ik_df[, c("shrub", "perennial", "annual", "herbaceous")],
                2, scale )
#check
head( ik_sc)

##########    save relevant data and workspaces     ###########
# Since we created multiple objects that we want to use in our analysis #
# we save the entire workspace and load it into the analysis script #
save.image( "MarkPrepWorkspace.RData" )
########## End of saving section ##################################
################## Save your data and workspace ###################

############# END OF SCRIPT ########################################