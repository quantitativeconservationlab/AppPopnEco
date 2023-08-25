#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                  ###
### Our study species are Great Horned Owl and Barn Owl in coastal  ###
### Texas. We used 3 call back surveys at x points spaced 1.6 km apart #
### to sample occupancy of these two species surrounding the area   ###
### where Aplomado Falcons were reintroduced, after becoming extinct # 
# in the U.S. in the 1940s. After multiple releases of captive-bred  #
## adult Aplomado Falcons in the 90s and early 2000s, recovery has been #
## stagnant. The working hypothesis is that encroachment of woody habitat #
## into open grasslands of Coastal Texas is bringing predatory owls into #
#  greater contact with Aplomado Falcons.                             ##
## Here we examine associations between owls' occupancy and habitat   ##
### for detection we look at the effects of wind speed, cloud cover,   ##
### and day of the year.                                               ##
###                                                                   ###
#########################################################################
########### clean workspace and load required packages ####################
###########################################################################

# Install new packages from "CRAN" repository. # 
install.packages( "jagsUI" ) #actually a collection of packages 

####### load relevant packages ###
library( dplyr ) #dataframe manipulations.
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( tidyr ) #to use spread and other df functions
library( ggplot2 ) #visualize data
library( jagsUI ) #to run RJAGS

########## end of package loading ###########
#######################    import relevant data   ##############
#####clean workspace to improve efficiency: ###
rm(list = ls() ) 
# Check that you are in the right project folder
getwd()

# set directory where your data are:
datadir <- paste( getwd(), "/Data/Texas/", sep = "" )

# load habitat predictors for the ecological model
habdf <- read.csv( file = paste( datadir, "habitat.csv",
              sep = ""), header = TRUE )

# load detection predictors
detdf <- read.csv( file = paste( datadir, 
              "TexasOwls_DetectionPredictors.csv", sep = ""),
                      header = TRUE )

#load Great Horned Owl detections
greatdf <- read.csv( file = paste( datadir, 
          "GreatHornedOwl_Detections.csv", sep = ""),
                   header = TRUE )

#load Barn Owl detections
barndf <- read.csv( file = paste( datadir, 
            "BarnOwl_Detections.csv", sep = ""),
                    header = TRUE )

####################################################################
##### Ready data for analysis --------------
### select which owl you want to analyse data for #
# Here we choose Great horned owls. For homework use Barn owls. 
#start by viewing detection dataframe
head( greatdf ); tail( greatdf); dim( greatdf )
# What does this tell you about the detection data?
#Answer: 
#

#create a dataframe to hold detections only 
y_obs <- greatdf[ ,c("p1", "p2", "p3") ]

#check 
head(y_obs)

#view habitat predictor data frame 
head( habdf )

#We explore predictors
#start by choosing which you want to look at
prednames <- c("grass", "forest", "shrub", "ah" ) 

# Check correlation among predictors.
round(cor( habdf[ , prednames] ),2)
# Are there any predictors we need to worry about?
# What correlations would be worrisome?
## Answer:
###

# loop over each to create histograms for each predictor:
for( p in 1:length(prednames) ){
  # create an object with the ggplot so that you can display it 
  # in a loop 
  a <- ggplot( habdf ) + #choose your data
    theme_bw( base_size = 15 ) + #choose a preset theme
    labs( x = prednames[p] ) + #label x axis using our predictor names
    geom_histogram( aes( get(prednames[p]) ), bins= 10 ) #plot histogram
  # display your plot object
  print( a )
}
# What do you note? Any apparent issues with these data?
# Answer:
#

#we are using habitat as static (no annual changes) 
#so we replicate values for both years by appending it to 
#the owl dataframe
habdf <- greatdf %>% 
  dplyr::select( Site.ID, yearid, siteyear_id ) %>% 
  left_join( habdf, by = "Site.ID" )
#view
head(habdf); dim(habdf)

#from these habitat predictors, we choose shrub cover (%), #
# and 'ah' which is aplomado habitat, a composite of open #
# habitats where we expect to find Aplomado Falcons. # 
# We extract those habitats into a new dataframe
X <- habdf[,c("shrub", "ah") ]
#standardise these predictors
X <- apply( X, MARGIN = 2, FUN = scale )
#view
head( X )

#Now explore detection predictors 
#view
tail(detdf);dim(detdf)
str(detdf)

#create a vector of each predictor for easy plotting
w <- as.vector( as.matrix(detdf[ ,c("wind1", "wind2", "wind3") ]))
c <- as.vector( as.matrix(detdf[ ,c("cloud1", "cloud2", "cloud3") ]))

#derive a histogram
hist( w )
hist( c )
#what does this tell you about which predictors we should use?
#what else could be relevant that we are ignoring?
#Answer:
# 

#we create a siteXsurvey matrix for each one and standardise them
#wind
wind_sc <- scale( detdf[ ,c("wind1", "wind2", "wind3") ])
#replace missing values with mean 0
wind_sc[is.na(wind_sc)] <- 0


#number of sites X year
I <- max( detdf$siteyear_id )
#number of replicate surveys each year
J <- 3


#########################################################################################
###########################################################################
####################### define MCMC settings ##############################

nt <- 10; nb <- 20000; nc <- 5 #thinning, burnin, chains

##### end of MCMC parameters definition #############
############################################################################
#################    run alternative models ################################
#############################################################################
#############################################################################
####### single-season single species occupancy model    ######
#############################################################################
############## Specify model in bugs language:  #####################
sink( "m1.txt" )
cat( "
     model{
     
      #priors
      #define intercept for occupancy
      # note that this is using the logit transformation #
      # do you recognise it?
      # it allows us to set the prior on the real scale #
      # where the intercept represents the mean probability averaged 
      # across sites and years...here we use a beta distribution 
      # which ranges between 0 to 1, the parameters given to this prior
      # should result in a tight, bell-shaped distribution centered around
      # 0.5 - 
      # if unsure you can plot different distributions to see what they look like
      
      #loop over years so that we get a year specific intercept
      for( k in 1:2 ){ 
        int.psi[ k ] <- log( mean.psi[ k ] / ( 1 - mean.psi[ k ] ) )
        mean.psi[ k ] ~ dbeta( 4, 4 )  #mean occupancy probability 
      }
        #We use the same process to define the prior
        # in the detection model
        int.p <- log( mean.p / ( 1 - mean.p ) )
        mean.p ~ dbeta( 4, 4 )  #mean detection prob
        
      #priors for fixed coefficients in ecological model:
      for( b in 1:B ){ #loop over number of predictors
      
        beta.psi[ b ] ~ dnorm( 0, 0.05 ) 
      
        #this precision is not fully uninformative but actually 
        #regularizes the estimates around 0 to help computation
      }
      
      #priors for fixed coefficients in detection submodel:
        alpha.p ~ dnorm( 0, 0.05 ) 

      #Define the ecological model for occupancy:
      for( i in 1:I ){  #loop over I number of sites
            
            #true occupancy state, z, is given a bernoulli distribution
            # with probability psi:
            z[ i ] ~ dbern( psi[ i ] ) 
            
            #probability of occupancy, psi is linked to predictors
            # using a logit function:
            logit( psi[ i ] ) <- int.psi[ yearid[i] ] + #year effect
                  beta.psi[ 1 ] * X[ i, 1 ] + #shrub
                  beta.psi[ 2 ] * X[ i, 2 ] #aplomado habitat
                  
        } #I sites
        
      #Define the observation model: 
      for( i in 1:I ){  #loop over sites
        for( j in 1:J ){ #loop over surveys each year
          
            # Link probability of detection, p, to predictors 
            # using a logit function
            logit( p[ i, j ] ) <- int.p +
                #fixed predictors
                alpha.p * wind_sc[ i, j ]  # wind
            
        #Here we link our observations to the estimated, true occupancy, z,
        # from our ecological model above
      
        y_obs[ i,j ] ~ dbern( z[ i ] * p[ i, j ] ) 

        #Estimate what the model would have produced as observations
        # we do this for model evaluation later
        yhat[ i,j ] ~ dbern( z[ i ] * p[ i, j ] )
        
        #Estimate the likelihood of observed and predicted
        # values for model validation later
        lik_yobs[ i,j ] <- ( ( psi[ i ] * p[ i,j ] )^y_obs[ i,j ] ) *
              ( ( 1 - psi[ i ] * p[ i,j ] )^( 1 - y_obs[ i,j ] ) )
        
        #likelihood of estimated detections:
        lik_yhat[ i,j ] <- ( ( psi[ i ]* p[ i,j ] )^yhat[ i,j ] ) *
            ( ( 1 - psi[ i ] * p[ i,j ] )^( 1 - yhat[ i,j ] ) )

      }}#
        
     } #model close
     
     ", fill = TRUE )

sink()
################ end of model specification  #####################################
modelname <- "m1.txt"
#parameters monitored #only keep those relevant for model comparisons
params <- c( 'int.psi' #intercept for occupancy model
             , 'int.p' #intercept for detection
             , 'beta.psi' #fixed coefficients for occupancy
             , 'alpha.p' #fixed coefficients for detection
             ,'z' #estimated occupancy state
             ,'psi' #probability of occupancy
             ,'p' #probability of detection
              , 'lik_yobs' #likelihood for each occupancy observation
              , 'lik_yhat' #likelihood for occupancy observations predicted by the model
            , 'yhat' #estimated occurrence from model
             
)

#create initial values for the model coefficients
zst <- rep(1, I )
#number of occupancy predictors
B <- 2
#number of detection predictors
A <- 1
#create initial values to start the algorithm
inits <- function(){ list( beta.psi = rnorm( B ),
                           alpha.p = rnorm( A ) 
                          ,z = zst
) }
#combine data into object:
str( great.data <- list( y_obs = y_obs, #observed occupancy for each species
                       J = J, I = I, A = A, B = B,
                       yearid = greatdf$yearid  
                       #site level predictors
                       ,X = as.matrix( X )
                       ,wind_sc = wind_sc
) )                

#call JAGS and summarize posteriors:
m_great <-  autojags( great.data, inits = inits, params, modelname, #
                 n.chains = nc, n.thin = nt, n.burnin = nb,
                 iter.increment = 20000, max.iter = 500000, 
                 Rhat.limit = 1.02,
                 save.all.iter = FALSE, parallel = TRUE ) 

###### end m1 ########
##################################################################
### save workspace ###
save.image( "OccBayesResults.RData" )

################# end of script #############################################