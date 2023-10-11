#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we import our cleaned data for the time-series of point count #
#  observations for Piute ground squirrels at the NCA and run        ##
## robust population N-mixture analyses. The models are hierarchical  #
#  with : (1) an ecological submodel linking abundance to             #
## environmental predictors; (2) an observation submodel linking     ##
##  detection probability to relevant predictors.                    ##
##                                                                   ##
# Female Piute ground squirrels give birth to an average of 5-10 young#
# Reproduction and survival are likely influenced by colder temperature #
# in Feb, when they come out of hibernation.                          #                
# Survival is likely affected by hot temperatures, with individuals   #
# unable to forage when temperatures are too hot.                     #
# Survival is expected to be higher in sites with more sagebrush.     #
#                                                                     #
# Detection may be related to observer effects and to time of day as a #
# quadratic, with higher detection expected in the middle of the day, #
# when squirrels are most active.                                     #
#                                                                    ##
#######################################################################
##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

# load packages:
library( tidyverse )#includes dplyr, tidyr and ggplot2
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( jagsUI )
## end of package load ###############
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Data/", sep = "" )
# load our cleaned data
rawdf <- read.csv( file = paste( datadir, "open_counts.csv", sep = ""),
                    header = TRUE )
#### End of data load -------------
####################################################################
##### Ready data for analysis --------------
#create new dataframe from imported data
opendf <- rawdf
#view
head( opendf ); dim( opendf ) 

# extract broad parameters of interest
#number of sites:
I <- length( unique(opendf$o.sites ))
#number of repeat visits
J <- 3
#number of primary surveys
K <- length(unique(opendf$year) )
#number of rows in the opendf
T <- dim(opendf)[1]

#we need site ids that are ordered starting from 1
opendf$siteid <- as.numeric(as.factor(opendf$o.sites))

#we need yearid also ordered starting from 1
opendf$yearid <- as.numeric(as.factor(opendf$year))

# extract site by year level predictors (including site and year id)
X <- opendf %>%  
  dplyr::select( siteid, yearid, cheatgrass, sagebrush, 
                 Feb.minT, AprMay.maxT )
#define predictor names
prednames <- c( "cheatgrass", "sagebrush", 
                "Feb.minT", "AprMay.maxT" )
#scale predictors
X[,prednames] <- apply( X[,prednames], MARGIN = 2, FUN = scale ) 

head( X )

#extract observation level predictors
#create index of columns of interest
tidx <- grep("time", colnames( opendf), value = FALSE)
oidx <- grep( "obs", colnames(opendf), value = FALSE)
#also for our observations
yidx <- grep( "count.j", colnames( opendf ), value = FALSE)

#now use those indices to automatically select correct columns and scale
time_sc <- scale( opendf[ ,tidx] )
#we want a quadratic term for time also
time2_sc <- scale( (opendf[ ,tidx])^2 )
# for observers, we will include them as a random intercept instead
# of a fixed effect so we convert them to a number to represent index
obvs <- as.matrix( opendf[,oidx] )
obvs <- as.numeric( as.factor(obvs) )
#turn back into matrix
obvs <- structure( obvs, dim = dim( opendf[,oidx] ), class = "matrix" )
#number of observers 
S <- max( obvs, na.rm = TRUE )

#extract maximum counts for year 1 for each site for initial lambda
maxcounts <- opendf %>% 
  #filter to year 1
  dplyr::filter( yearid == 1 ) %>% 
  #select count columns
 dplyr::select( count.j1, count.j2, count.j3 ) 
head( maxcounts )
#estimate max for each site
maxcounts <- apply( maxcounts, MARGIN = 1, FUN = max )

#It really helps to provide initial values for N. Here we make those
#from our max counts at each site (what traditionally would have been 
# the index of abundance)
Ninits <- opendf %>% 
  #select count columns
  mutate( maxc  = count.j1 + count.j2 + count.j3 ) %>% 
  dplyr::select( siteid, yearid, maxc ) 

head(Ninits)
### cannot have 0 abundance so we replace those with 1
Ninits$maxc[ which(Ninits$maxc == 0)] <-1
Ninits <- pivot_wider( Ninits, names_from = yearid, values_from = maxc)
Ninits <- as.data.frame( Ninits ) %>% 
  select( -siteid )
Ninits
dim(Ninits)

##########################################################################
####### m1 N-mixture abundance model assuming N ~ Poisson(gamma)   #
# where gamma takes the traditional role of lambda                  #
# gamma predictors: cheatgrass, sagebrush, Feb.minT, AprMay.maxT #
# detection predictors: time, time^2 and observer as random intercept #
############################################################################
############## Specify model in bugs language:  #####################
sink( "m1.txt" )
cat( "
     model{
     
      #priors
      #for detection model: 
      #define intercept as mean probs:
      int.det <- log( mean.det / ( 1 - mean.det ) )
      mean.det ~ dbeta( 4, 4 )
      
      #random intercept for observer
      for( s in 1:S ){
        eps.det[s] ~ dnorm( 0, pres.det ) T(-7, 7)
      }
      #associated variance of random intercepts:     
      pres.det <- 1/ ( sigma.det * sigma.det )
      #sigma prior specified as a student t half-normal:
      sigma.det ~ dt( 0, 2.5, 7 ) T( 0, )
      
      #priors for detection coefficients:
      #define as a slightly informative prior
      for( a in 1:A){
        alpha[ a ] ~ dnorm( 0, 0.2 ) T(-7, 7 )
      }
      #priors for abundance coefficients:
      for( b in 1:B ){
        #define as a slightly informative prior
        beta[ b ] ~ dnorm( 0, 0.2 ) T(-7, 7 )
      }
      
      #prior for abundance model intercept 
      int.gam ~ dgamma( 0.01, 0.01 )

     # ecological model of abundance
    for( i in 1:I ){

      #Abundance in year 1 is just a Poisson with lambda derived
      #from maximum counts observed at that site
       N[ i,1 ] ~ dpois( maxcounts[i] )

       for( k in 2:K ){
       #Abundance in following years is related to gamma
       N[ i, k ] ~ dpois(  gamma[ i, k ] )
    
        } #close K
     } #close I
    
    #observation model
      for( t in 1:T ){ #loop over each row of your opendf
      
         #link gamma to predictors
          log( gamma[ siteid[t], yearid[t] ] ) <- int.gam +
                                    inprod( beta, X[ t, ]  )
         
        for( j in 1:J ){ #loop over surveys

        #model probability of detection
        logit( p[siteid[t], yearid[t], j ] ) <- int.det + 
                  #random intercept for observer effect
                  eps.det[ obvs[t,j] ] +
                  #quadratic effect of time of day
                  alpha[1] * time[ t, j ] +
                  alpha[2] * time2[ t, j ]
                  
        #observed counts distributed as a Binomial:
        y_obs[ t, j ] ~ dbin( p[ siteid[t], yearid[t], j ], 
                              N[ siteid[t], yearid[t] ] )  
                  
    } #close J
    } #close T
    
    } #model close
     
     ", fill = TRUE )

sink()
################ end of model specification  #####################################     
modelname <- "m1.txt"
#parameters monitored
params <- c(  'int.det' #intercept for detection
              , 'int.gam' #intercept for lamda
              , 'alpha' #detection coefficients
              , 'eps.det' #random intercepts in detection
              , 'sigma.det' #error for random intercept
              , 'beta' #abundance coefficients
              , 'gamma' #abundance rate 
              , 'p' #estimate of detection probability
              , 'N' #estimates of abundance

)

#how many ecological predictors that are fixed effects
B <- 4
#how many detection predictors that are fixed effects
A <- 2
#define initial parameter values
inits <- function(){ list( beta = rnorm( B ),
                           alpha = rnorm( A ),
                           N = as.matrix(Ninits) ) }

#define data that will go in the model
str( win.data <- list( y_obs = as.matrix( opendf[ ,yidx] ),
                       #number of sites, surveys, det predictors, and abund preds
                       I = I, J = J, A = A, B = B, S = S, K= K, T = T,
                       siteid = opendf$siteid,
                       yearid = opendf$yearid,
                       #max counts
                       maxcounts = maxcounts,
                      #siteXyear level habitat predictors
                       X = X[,c("cheatgrass", "sagebrush", "Feb.minT", "AprMay.maxT")],
                       #observation predictors:
                       time = time_sc,
                       time2 = time2_sc,
                       obvs = obvs
) )

#call JAGS and summarize posteriors:
m1 <- autojags( win.data, inits = inits, params, modelname, #
                n.chains = 5, n.thin = 10, n.burnin = 0,
                iter.increment = 10000, max.iter = 500000, 
                Rhat.limit = 1.05,
                save.all.iter = FALSE, parallel = TRUE ) 

thanks#view results 
summary(m1)
plot(m1)
############################################################################
##########################################################################
####### m2 Gompertz N-mixture abundance model    #
# We extend model 1 to incorporate a gompert process that allows for #
# density dependence and estimates of instantaneous population growth #
# rate. 
# ecological predictors: cheatgrass, sagebrush, Feb.minT, AprMay.maxT #
# detection predictors: time, time^2 and observer as random intercept #
############################################################################
############## Specify model in bugs language:  #####################
sink( "m2.txt" )
cat( "
     model{
     
      #priors
      #for detection model: 
      #define intercept as mean probs:
      int.det <- log( mean.det / ( 1 - mean.det ) )
      mean.det ~ dbeta( 4, 4 )
      
      #random intercept for observer
      for( s in 1:S ){
        eps.det[s] ~ dnorm( 0, pres.det ) T(-7, 7)
      }
      #associated variance of random intercepts:     
      pres.det <- 1/ ( sigma.det * sigma.det )
      #sigma prior specified as a student t half-normal:
      sigma.det ~ dt( 0, 2.5, 7 ) T( 0, )
      
      #priors for detection coefficients:
      #define as a slightly informative prior
      for( a in 1:A){
        alpha[ a ] ~ dnorm( 0, 0.2 ) T(-7, 7 )
      }
      #priors for abundance coefficients:
      for( b in 1:B ){
        #define as a slightly informative prior
        beta[ b ] ~ dnorm( 0, 0.2 ) T(-7, 7 )
      }
      
      #prior for abundance model intercept 
      int.gam ~ dgamma( 0.01, 0.01 )
      # int.gam ~ dnorm( log(1.3), 0.5 )
      
      #prior for omega (equilibrium abundance)
      omega ~ dunif( 40, 100 )
      
    # ecological model of abundance
    for( i in 1:I ){

      #Abundance in year 1 is just a Poisson with lambda derived
      #from maximum counts observed at that site
       N[ i,1 ] ~ dpois( maxcounts[i] )

       for( k in 2:K ){

       #define abundance from a Gompertz autoregressive process
       # where gamma is instantaneous population growth rate
       # and omega is equilibrium abundance (K)
          lambda[ i, k ] <- N[ i, k-1 ] * gamma[ i, k ] *
                ( 1 - log( N[ i, k-1 ] ) / log( omega + 1 ) ) 
          N[ i, k ] ~ dpois( lambda[ i, k] )
          
        } #close K
     } #close I
    
    #observation model
      for( t in 1:T ){ #loop over each row of your opendf
      
         #link gamma to predictors
          log( gamma[ siteid[t], yearid[t] ] ) <- int.gam +
                                    inprod( beta, X[ t, ]  )
         
        for( j in 1:J ){ #loop over surveys

        #model probability of detection
        logit( p[siteid[t], yearid[t], j ] ) <- int.det + 
                  #random intercept for observer effect
                  eps.det[ obvs[t,j] ] +
                  #quadratic effect of time of day
                  alpha[1] * time[ t, j ] +
                  alpha[2] * time2[ t, j ]
                  
        #observed counts distributed as a Binomial:
        y_obs[ t, j ] ~ dbin( p[ siteid[t], yearid[t], j ], 
                              N[ siteid[t], yearid[t] ] )  
                  
    } #close J
    } #close T
    
    } #model close
     
     ", fill = TRUE )

sink()
################ end of model specification  #####################################     
modelname <- "m2.txt"
#parameters monitored
params <- c(  'int.det' #intercept for detection
              , 'int.gam' #intercept for lamda
              , 'alpha' #detection coefficients
              , 'eps.det' #random intercepts in detection
              , 'sigma.det' #error for random intercept
              , 'beta' #abundance coefficients
              , 'gamma' #maximum instantaneous population growth rate (r) 
              , 'omega' #equilibrium abundance
              , 'p' #estimate of detection probability
              , 'N' #estimates of abundance
              
)

#how many ecological predictors that are fixed effects
B <- 4
#how many detection predictors that are fixed effects
A <- 2
#define initial parameter values
inits <- function(){ list( beta = rnorm( B ),
                           alpha = rnorm( A ),
                           N = as.matrix(Ninits) ) }

#define data that will go in the model
str( win.data <- list( y_obs = as.matrix( opendf[ ,yidx] ),
                       #number of sites, surveys, det predictors, and abund preds
                       I = I, J = J, A = A, B = B, S = S, K= K, T = T,
                       siteid = opendf$siteid,
                       yearid = opendf$yearid,
                       #max counts
                       maxcounts = maxcounts,
                       #siteXyear level habitat predictors
                       X = X[,c("cheatgrass", "sagebrush", "Feb.minT", "AprMay.maxT")],
                       #observation predictors:
                       time = time_sc,
                       time2 = time2_sc,
                       obvs = obvs
) )

#call JAGS and summarize posteriors:
m2 <- autojags( win.data, inits = inits, params, modelname, #
                n.chains = 5, n.thin = 10, n.burnin = 0,
                iter.increment = 10000, max.iter = 500000, 
                Rhat.limit = 1.05,
                save.all.iter = FALSE, parallel = TRUE ) 

#view results 
summary(m2)
plot(m2)
############################################################################

################## Save your data and workspace ###################
# Save workspace:
save.image( "RobCountBayesResults.RData" )

######### End of saving section ##################################

############# END OF SCRIPT #####################################

