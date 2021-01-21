#######################################################################
#######################################################################
####   This script was created by Dr. Jen Cruz as part of      ########
####          the Applied Population Ecology Class              #######
####                                                           ########
#### Their abundance is influenced by drought, low temperatures when ##
#### they emerge from hibernation in Feb and high temperatures in April-May #



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
