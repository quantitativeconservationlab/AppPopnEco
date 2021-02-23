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
# three days. Surveys were repeated over multiple years () #
# Surveys involved point counts for 2min where all individuals detected #
# over a x radius were recorded. 
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

# detection for point counts is related to time of day as a quadratic #

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
