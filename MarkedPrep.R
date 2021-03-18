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
#save.image( "CountPrepWorkspace.RData" )
########## End of saving section ##################################
################## Save your data and workspace ###################

############# END OF SCRIPT ########################################