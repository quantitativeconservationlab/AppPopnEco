#######################################################################
#######################################################################
####   This script was created by Dr. Jen Cruz as part of      ########
####          the Applied Population Ecology Class              #######
####                                                           ########
#### In this script we simulate the data that we will use      ########
#### for all the sample scripts provided in this course.       #######
#### We simulate three datasets: occupancy records for multiple ######
#### sites over multiple years. This will involve repeat surveys #####
###  within a primary season, when the populations are assumed to #####
###  be closed, with the polulations assumed to be open between  #####
#### primary seasons. This is the classic ROBUST DESIGN.         #####
### We will use a robust design to also simulate sampling over a #####
### subset of these sites using repeated counts. Lastly, another ####
### smaller subset of sites will be simulated to be sampled using ####
### capture-recapture methods, again following a robust design.   #####
####                                                             ##### 
### We link these demographic parameters (occupancy, abundance,   ###
### site persistence and extinction, survival ) to predictors    ####
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

# load packages:
library( tidyverse ) 
library( unmarked )

###################################################################
#### Load or create data -----------------------------------------
#### Simulating occupancy data #################
# Our study species is the Piute ground squirrel, Spermophilus    # 
# mollis. Ground squirrels are widely distributed in sagebrush-steepe #
# habitats of the Great Basin and Columbia Plateau. #
# Their abundance is influenced by drought, low temperatures when #
# they emerge from hibernation in Feb high temperatures in April-May #

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

# Number of sites we visit to search for their presence:
Io <- 100
# Number of sites we visit to count them without marking them:
Ic <- 50
# Number of sites we visit to trap them
Im <- 20


#############end of section creating data #########################
################## Save your data and workspace ###################

########## End of saving section ##################################
###################   END OF SCRIPT  ################################